//Convert vdif format to psrfits search mode format in 32-bit float
//Input vdif configuration is 2-bit real-sampled, 2048 MHz bandwdith, one frequency channel, two pols in separated files
//Sampling interval 1/2048/2 microsecond
//One frame 8192 bytes, thus 8 microsecond per frame (one frame contains only one pol)


#define VDIF_HEADER_BYTES       32

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <malloc.h>
#include "vdifio.h"
#include "vdif2psrfits.h"
#include "psrfits.h"
#include "dec2hms.h"
#include <mpi.h>
#include <fftw3.h>
#include <stdbool.h>

static uint32_t VDIF_BW = 2048; //Bandwidth in MHz
static uint32_t VDIF_BIT = 2;   //Bit per sample
static uint32_t VDIF_NCHAN = 1; //Number of channels

void usage(char *prg_name)
{
  fprintf(stdout,
	  "%s [options]\n"
	  " -f   Observing central frequency (MHz)\n"
          " -i   Input vdif pol0\n"
	  " -j   Input vdif pol1\n"
	  " -T   Telescope name (by default Pico Veleta)\n"
	  " -b   Band sense (-1 for lower-side, 1 for upper-side, by default 1)\n"
	  " -s   Seconds to get statistics to fill in invalid frames\n"
	  " -t   Time sample scrunch factor (by default 1). One time sample 8 microsecond\n"
	  " -S   Name of the source (by default J0835-4510)\n"
          " -r   RA of the source (AA:BB:CC.DD)\n"
          " -c   Dec of the source (+AA:BB:CC.DD)\n"
	  " -D   Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, S for Stokes, P for polarised signal, S for stokes, by default C)\n"
	  " -n   Number of channels kept (Power of 2 up to 4096, by default 1)\n"
	  " -v   Verbose\n"
	  " -O   Route of the output file \n"
	  " -h   Available options\n"
	  " Open MPI options\n"
	  " -p   Number of processors (by default from MPI_Comm_size)",
	  prg_name);
  exit(0);
}

// get offset number of VDIF frame from beginning of observation
int64_t getVDIFFrameOffset(const vdif_header *headerst, const vdif_header *header, uint32_t fps)
{
  uint32_t day[2],sec[2],num[2];
  int64_t offset;

  // Get epoch                                                                                 
  day[0]=getVDIFFrameMJD(headerst);
  day[1]=getVDIFFrameMJD(header);

  // Get second after epoch 
  sec[0]=getVDIFFrameSecond(headerst);
  sec[1]=getVDIFFrameSecond(header);

  // Get number of frame after second   
  num[0]=getVDIFFrameNumber(headerst);
  num[1]=getVDIFFrameNumber(header);

  offset=(day[1]-day[0])*86400*fps+(sec[1]-sec[0])*fps+(num[1]-num[0]);

  return offset;
}


int main(int argc, char *argv[])
{
  FILE *vdif[2],*out;
  bool pval[2], ifverbose, ifpol[2], ifout, chkend[2], pend[2];
  struct psrfits pf;
  
  char vname[2][1024],oroute[1024],ut[30],dat,vfhdr[2][VDIF_HEADER_BYTES],vfhdrst[VDIF_HEADER_BYTES],srcname[16],dstat,ra[64],dec[64],telname[64];
  int arg,n_f,i,j,k,fbytes,vd[2],nf_stat,ct,tsf,dati,nchan,npol,bs,Nts,nthd,nread[2],nf_skip;
  float freq,s_stat,fmean[2][2], *in_p0, *in_p1,s_skip, ***detmap;
  double mjd[2];
  long int idx[2], seed, chunksize[2],chunksize_org,Nfm,nfm_p[2],nskip, index[2];
  unsigned char *buffer[2], *obuffer[2], *chunk[2], *buffall[2];
  time_t t;
  double mean[4],sq,rms[j],spf;
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  int64_t offset_pre[2],offset[2];
  uint32_t fps,inval,inval_sub;

  // Set default values
  freq=0.0;
  s_stat=0.0;
  tsf=1;
  bs=1;
  strcpy(srcname,"J0835-4510");
  strcpy(telname,"Pico Veleta");
  dstat='C';
  npol=4;
  nchan=1;
  chunksize_org=1000000000;
  nthd=1;
  inval=0;
  ifverbose = false;
  ifout = false;
  for(i=0;i<2;i++) {
    ifpol[i] = false;
    chkend[i]=false;
    pend[i]=false;
    chunksize[i] = chunksize_org;
  }

  //Read arguments
  while ((arg=getopt(argc,argv,"hf:i:j:T:b:s:t:O:S:D:n:r:c:v")) != -1)
    {
      switch(arg)
	{
	case 'f':
	  freq=atof(optarg);
	  break;
		  
	case 'i':
	  strcpy(vname[0],optarg);
	  ifpol[0]=true;
	  break;

	case 'b':
	  bs=atoi(optarg);
	  break;

	case 'j':
	  strcpy(vname[1],optarg);
	  ifpol[1]=true;
	  break;

	case 's':
	  s_stat=atof(optarg);
	  break;
		  
        case 't':
	  tsf=atoi(optarg);
	  break;

	case 'r':
	  strcpy(ra,optarg);
	  break;

	case 'c':
	  strcpy(dec,optarg);
	  break;

	case 'T':
	  strcpy(telname,optarg);
	  break;

	case 'S':
	  strcpy(srcname,optarg);
	  break;
		  
        case 'D':
	  strcpy(&dstat,optarg);
	  if(dstat!='C') npol=1;
	  break;
		  
	case 'O':
	  strcpy(oroute,optarg);
	  ifout=true;
	  break;

	case 'n':
	  nchan=atoi(optarg);
	  break;

	case 'v':
	  ifverbose=true;
	  break;
		  
	case 'h':
	  usage(argv[0]);
	  return 0;

	default:
	  usage(argv[0]);
	  return 0;
	}
    }
  
  //Check if arguments are enough to procceed
  if(freq==0.0)
	{
	  fprintf(stderr,"Missing info of observing central frequency.\n");
	  exit(0);
	}
  
  if(ifpol[0] == false)
	{
	  fprintf(stderr,"No input file provided for pol0.\n");
	  exit(0);
	}
  
  if(ifpol[1] == false)
	{
	  fprintf(stderr,"No input file provided for pol1.\n");
	  exit(0);
	}

  if(tsf<1)
	{
	  fprintf(stderr,"Invalid sample scrunching factor.\n");
	  exit(0);
	}  
  
  if(ifout == false)
	{
	  fprintf(stderr,"No output route specified.\n");
	  exit(0);
	}

  float det[nchan][4],sdet[nchan][4];
  
  //Get seed for random generator
  srand((unsigned)time(&t));
  seed=0-t;

  //Read the first header of vdif pol0
  vdif[0]=fopen(vname[0],"rb");
  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
  fclose(vdif[0]);
  
  //Get header info
  /*-----------------------*/
  //Get frame bytes
  fbytes=getVDIFFrameBytes((const vdif_header *)vfhdr[0])-VDIF_HEADER_BYTES;
  printf("Bytes per frame: %i\n",fbytes);

  //Calculate time interval (in unit of microsecond) of a frame
  spf=(double)fbytes/VDIF_BIT*8/(VDIF_BW*2);

  //Frame per second
  fps=1000000*VDIF_BIT/8*2*VDIF_BW/fbytes;

  //Calculate how many frames to get statistics
  nf_stat=s_stat*1.0e6/spf;

  //Calculate number of frames in one read chunk
  Nfm=chunksize_org/(fbytes+VDIF_HEADER_BYTES);
  printf("Number of frame per reading chunk: %ld\n",Nfm);

  // Calculate how many frames to skip from the beginning
  nf_skip=s_skip*1.0e6/spf;
  printf("Number of frames to skip from the beginning: %i.\n",nf_skip);

  //Allocate memo for frames
  for(j=0;j<2;j++)
	{
	  buffer[j]=(unsigned char *)malloc(sizeof(unsigned char)*fbytes);
	  obuffer[j]=(unsigned char *)malloc(sizeof(unsigned char)*fbytes*4);
	}

  // Prepare FFT
  Nts = fbytes*4;
  in_p0 = (float *) malloc(sizeof(float)*Nts);
  in_p1 = (float *) malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  printf("Determining FFT plan...length %d...",Nts);
  pl0 = fftwf_plan_dft_r2c_1d(Nts, in_p0, out_p0, FFTW_MEASURE);
  pl1 = fftwf_plan_dft_r2c_1d(Nts, in_p1, out_p1, FFTW_MEASURE);
  if (pl0 != NULL && pl1 != NULL) 
    printf("Done.\n");
  else
    {
      fprintf(stderr,"Error in creating FFT plan.\n");
      exit(0);
    }

  // Scan the beginning specified length of data, choose valid frames to get mean of total value in each frame
  fprintf(stderr,"Scan %.2f s data to get statistics...\n",s_stat); 
  for(j=0;j<2;j++)
    {
      mean[j]=0.0;
      sq=0.0;
      rms[j]=0.0;
      ct=0;

      // Open file
      vdif[j]=fopen(vname[j],"rb");

      fseek(vdif[j],(VDIF_HEADER_BYTES+fbytes)*nf_skip,SEEK_CUR);

      for(i=0;i<nf_stat;i++)
		{
		  //Read header
		  fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);
		  
		  //Valid frame
		  if(!getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[j],VDIF_HEADER_BYTES+fbytes,ifverbose))
			{
			  //Read data in frame
			  fread(buffer[j],1,fbytes,vdif[j]);

			  //Expand from 2-bit to 8-bit
			  convert2to8(obuffer[j], buffer[j],fbytes);

			  //Accumulate values
			  for(k=0;k<fbytes*4;k++)
				{
				  dati=(int)obuffer[j][k];
				  mean[j]+=(double)dati;
				  sq+=pow((double)dati,2.0);
				}

			  //Add counter
			  ct++;
			}
		  //Invalid frame
		  else
			{
			  //Skip the data
			  fseek(vdif[j],fbytes,SEEK_CUR);
			}
		}
	  fclose(vdif[j]);
	  mean[j]=mean[j]/ct/fbytes/4;
	  rms[j]=sqrt(sq/ct/fbytes/4-pow(mean[j],2.0));
	  printf("Pol%i: mean %lf, rms %lf.\n",j,mean[j],rms[j]);
	}

  // Open VDIF files for data reading
  for(j=0;j<2;j++)
    vdif[j]=fopen(vname[j],"rb");

  // Calibrate the difference in starting time between the two pols
  printf("Calibrating potential difference in starting time between two pols...\n");

  // Move to the first valid frame for both pols
  for(j=0;j<2;j++)
    {
      // Move to the first valid frame and get header info
      do {
	// Get header
	fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);

	// Valid frame
	if(!getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[j],VDIF_HEADER_BYTES+fbytes, ifverbose))
	  {
	    mjd[j]=getVDIFFrameDMJD((const vdif_header *)vfhdr[j], fps);
	    break;
	  }
	  // Invalid frame
	else
	  fseek(vdif[j],fbytes,SEEK_CUR);
      }while(feof(vdif[j])!=1);
    }

  // Synchronize starting time
  while(mjd[0]!=mjd[1] || getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[0],VDIF_HEADER_BYTES+fbytes,ifverbose) || getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[1],VDIF_HEADER_BYTES+fbytes,ifverbose))
    {
      if(mjd[0]<mjd[1])
	j=0;
      else
	j=1;

      fseek(vdif[j],fbytes,SEEK_CUR);
      fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);
      mjd[j]=getVDIFFrameDMJD((const vdif_header *)vfhdr[j], fps);
    }

  //Get starting MJD and UT
  memcpy(vfhdrst,vfhdr[0],VDIF_HEADER_BYTES);
  mjd2date(mjd[0],ut);
  printf("Starting time synchronized. Start UT of VDIF: %s\n",ut);

  // Set starting position for reading
  for(j=0;j<2;j++)
    {
      fseek(vdif[j],-VDIF_HEADER_BYTES,SEEK_CUR);
      offset_pre[j]=-1;
    }

  //Set psrfits main header
  printf("Setting up PSRFITS output...\n");
  pf.filenum = 0;           // This is the crucial one to set to initialize things
  pf.rows_per_file = 200;  // Need to set this based on PSRFITS_MAXFILELEN

  //Set values for our hdrinfo structure
  pf.hdr.scanlen = 86400; // in sec
  strcpy(pf.hdr.observer, "A. Eintein");
  strcpy(pf.hdr.telescope, telname);
  strcpy(pf.hdr.obs_mode, "SEARCH");
  strcpy(pf.hdr.backend, "MARK6");
  strcpy(pf.hdr.source, srcname);
  strcpy(pf.hdr.date_obs, ut);
  strcpy(pf.hdr.poln_type, "LIN");
  // Specify status of output data
  if( dstat == 'C' )
    strcpy(pf.hdr.poln_order, "AABBCRCI");
  else if ( dstat == 'S' )
    strcpy(pf.hdr.poln_order, "IQUV");
  else if ( dstat == 'P' )
    strcpy(pf.hdr.poln_order, "AABB");
  else
    strcpy(pf.hdr.poln_order, "AA+BB");
  strcpy(pf.hdr.track_mode, "TRACK");
  strcpy(pf.hdr.cal_mode, "OFF");
  strcpy(pf.hdr.feed_mode, "FA");
  pf.hdr.dt = spf*tsf/1.0e6;
  printf("Sampling interval: %0.9lf\n",pf.hdr.dt);
  // Shift central frequency half a raw spectral width up from real fft
  pf.hdr.fctr = freq+1.0/spf/2;
  pf.hdr.BW = VDIF_BW;
  pf.hdr.nchan = nchan;
  pf.hdr.MJD_epoch = mjd[0];
  strcpy(pf.hdr.ra_str,ra);
  strcpy(pf.hdr.dec_str,dec);
  pf.hdr.azimuth = 123.123;
  pf.hdr.zenith_ang = 23.0;
  pf.hdr.beam_FWHM = 0.25;
  pf.hdr.ibeam = 0;
  pf.hdr.start_lst = 10000.0;
  pf.hdr.start_sec = 25000.82736876;
  pf.hdr.start_day = 55000;
  pf.hdr.scan_number = 1;
  pf.hdr.rcvr_polns = 2;
  pf.hdr.onlyI = 0;
  pf.hdr.summed_polns = 0;
  pf.hdr.offset_subint = 0;
  pf.hdr.orig_nchan = pf.hdr.nchan;
  pf.hdr.orig_df = pf.hdr.df = pf.hdr.BW / pf.hdr.nchan;
  pf.hdr.nbits = 32;
  pf.hdr.npol = npol;
  pf.hdr.chan_dm = 0.0;
  pf.hdr.fd_hand = 1;
  pf.hdr.fd_sang = 0;
  pf.hdr.fd_xyph = 0;
  pf.hdr.be_phase = 1;
  pf.hdr.nsblk = 12500;
  pf.hdr.ds_time_fact = 1;
  pf.hdr.ds_freq_fact = 1;
  sprintf(pf.basefilename, "%s/%s",oroute,ut);

  psrfits_create(&pf);

  //Set values for our subint structure
  pf.sub.tsubint = pf.hdr.nsblk * pf.hdr.dt;
  pf.tot_rows = 0;
  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
  pf.sub.lst = pf.hdr.start_lst;
  pf.sub.ra = pf.hdr.ra2000;
  pf.sub.dec = pf.hdr.dec2000;
  pf.sub.feed_ang = 0.0;
  pf.sub.pos_ang = 0.0;
  pf.sub.par_ang = 0.0;
  pf.sub.tel_az = pf.hdr.azimuth;
  pf.sub.tel_zen = pf.hdr.zenith_ang;
  pf.sub.bytes_per_subint = (pf.hdr.nbits * pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk) / 8;
  pf.sub.FITS_typecode = TBYTE;  // 11 = byte

  // Create and initialize the subint arrays
  pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  for (i = 0 ; i < pf.hdr.nchan ; i++)
	{
	  pf.sub.dat_freqs[i] = pf.hdr.fctr - 0.5 * pf.hdr.BW + 0.5 * pf.hdr.df + i * pf.hdr.df;
	  pf.sub.dat_weights[i] = 1.0;
	}
  pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pf.sub.dat_scales = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  for (i = 0 ; i < pf.hdr.nchan * pf.hdr.npol ; i++)
	{
	  pf.sub.dat_offsets[i] = 0.0;
	  pf.sub.dat_scales[i] = 1.0;
	}

  pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);

  printf("Header prepared. Start to write data...\n");

  // First read of data chunk
  for(i=0;i<2;i++)
    {
      chunk[i] = (unsigned char *)malloc(chunksize[i]+fbytes+VDIF_HEADER_BYTES);
      nread[i] = fread(chunk[i],1,chunksize[i],vdif[i]);
      nfm_p[i] = Nfm;
      index[i] = 0;
    }

  // Allocate memory for detection map and filler judge
  // Get 2 times the array size to deal with potential fake samples
  nts_map = chunksize_org / (fbytes + VDIF_HEADER_BYTES) * 2;
  float detmap[nts_map][nchan][4];
  int valtsmarker[nts_map];
  for(i=0;i<nts_map;i++)
    for(j=0;j<nchan;j++)
      for(k=0;k<4;k++)
	{
	  detmap[i][j][k] = 0.0;
	}

  // MPI environment setup
  MPI_Init(NULL, NULL);

  int nthreads, MPIi, world_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  unsigned char *recvbuff[0] = malloc(sizeof(unsigned char) * detsum * fbytes);
  unsigned char *recvbuff[1] = malloc(sizeof(unsigned char) * detsum * fbytes);
  unsigned char *MPIbuff[0] = malloc(sizeof(unsigned char) * fbytes);
  unsigned char *MPIbuff[1] = malloc(sizeof(unsigned char) * fbytes);
  float MPIdet[nchan][4], *detarr;

  if(world_rank == 0)
    {
      // Sum of detections per thread
      detsum = chunksize_org / (fbytes + VDIF_HEADER_BYTES) / nthreads;

      // Buffer for all valid detections in one loop
      for(i=0;i<2;i++)
	buffall[i] = (unsigned char *)malloc(fbytes*detsum*nthreads);

      // Return buffer of detections
      detarrglobal = (float *)malloc(sizeof(float)*detsum*nchan*4*nthreads);
    }

  // Main loop to write subints
  do
    {
      /* Prepare data buffer for detection jobs on head node */
      if(world_rank == 0)
	{
	  inval_sub = 0;
	  tscount = 0;
	  valtscount = 0;

	  // Store buffer for each valid detection
	  while(valtscount < detsum * nthreads)
	    {
	      // Consecutive check on both pols
	      for(j=0;j<2;j++)
		{
		  // Get real frame header while dealing with non-integer gaps in between frames
		  nskip=0;
		  for(;;)
		    {
		      // Not enough data in chunk to fill header
		      if(index[j] + VDIF_HEADER_BYTES > chunksize[j])
			{
			  if(ifverbose)
			    fprintf(stdout,"Pol %i read new chunk.\n",j);

			  // Move leftover to the beginning
			  memmove(chunk[j],chunk[j]+index[j],chunksize[j]-index[j]);

			  //  If not end of pol file, read another chunk
			  if(!chkend[j])
			    nread[j]=fread(chunk[j]+chunksize[j]-index[j],1,chunksize_org, vdif[j]);
			  // End of pol file, break
			  else
			    {
			      pend[j]=true;
			      break;
			    }

			  // Mark it is the last read when at the end of pole file after read
			  if(nread[j]<chunksize_org)
			    chkend[j]=true;

			  // Set new chunk size and read position (index)
			  chunksize[j] = nread[j] + chunksize[j]-index[j];
			  index[j]=0;;
			}

		      // Shift forward by bytes until having a meaningful framebyte number
		      memcpy(vfhdr[j],chunk[j]+index[j],VDIF_HEADER_BYTES);
		      if (getVDIFFrameBytes((const vdif_header *)vfhdr[j]) == fbytes+VDIF_HEADER_BYTES)
			break;
		      else
			{
			  index[j]++;
			  nskip++;
			}
		    }

		  // Break when the end of pol file
		  if( pend[j] )
		    break;
	    
		  if(ifverbose && nskip>0)
		    fprintf(stdout,"Pol%i Skipped %i bytes.\n",j,nskip);

		  pval[j] = true;

		  // Get frame time offset
		  offset[j]=getVDIFFrameOffset((const vdif_header *)vfhdrst, (const vdif_header *)vfhdr[j], fps);

		  // Valid frame and Gap from the last frames
		  if( !getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[j],VDIF_HEADER_BYTES+fbytes,ifverbose) && offset[j] > offset_pre[j]+1 )
		    {
		      if(ifverbose)
			fprintf(stderr,"Pol%i: Current frame (%Ld) not consecutive from previous (%Ld).\n",j,offset[j],offset_pre[j]);
		      pval[j] = false;
		      fseek(vdif[j],-VDIF_HEADER_BYTES,SEEK_CUR);
		      offset_pre[j]++;
		    }
		  // Consecutive
		  else
		    {
		      // Not enough data in chunk to fill in frame
		      if(index[j] + VDIF_HEADER_BYTES + fbytes > chunksize[j])
			{
			  if(ifverbose)
			    fprintf(stdout,"Pol %i read new chunk.\n",j);

			  // Move leftover to the beginning
			  memmove(chunk[j],chunk[j]+index[j],chunksize[j]-index[j]);

			  // If not at end of pol file
			  if(!chkend[j])
			    nread[j]=fread(chunk[j]+chunksize[j]-index[j],1,chunksize_org, vdif[j]);
			  else
			    {
			      pend[j]=true;
			      break;
			    }

			  // Mark it is the last read
			  if(nread[j]<chunksize_org)
			    chkend[j]=true;

			  // New data chunk size and index 
			  chunksize[j] = nread[j] + chunksize[j]-index[j];
			  index[j]=0;
			}

		      // Get data in frame into buffer
		      memcpy(buffer[j], chunk[j] + index[j] + VDIF_HEADER_BYTES, fbytes);
		      offset_pre[j]++;
		      index[j]+=fbytes+VDIF_HEADER_BYTES;
		    }
		}

	      // Both pol consecutive
	      if(pval[0] == true && pval[1] == true)
		{
		  // Valid frame
		  if(!getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[0],VDIF_HEADER_BYTES+fbytes,ifverbose) && !getVDIFFrameInvalid_robust((const vdif_header *)vfhdr[1],VDIF_HEADER_BYTES+fbytes, ifverbose))
		    {
		      // Store data of the valid frame, leave marker
		      memcpy(buffall[0]+valtscount*fbytes, buffer[0], fbytes);
		      memcpy(buffall[1]+valtscount*fbytes, buffer[1], fbytes);
		      valtsmarker[valtscount] = tscount;
		      tscount++; valtscount++;
		    }
		  // Invalid frame
		  else
		    {
		      // Create fake detection with measured mean
		      if(ifverbose)
			fprintf(stderr,"Invalid frame detected in file %d subint %d (%f sec). Fake detection with measured mean.\n", pf.filenum, pf.tot_rows, pf.T);
		      getVDIFFrameFakeDetection_1chan(mean,rms,nchan,det,seed,fbytes,dstat);

		      for(j=0;j<nchan;j++)
			{
			  detmap[tscount][j][0] = dstat[j][0];
			  detmap[tscount][j][1] = dstat[j][1];
			  detmap[tscount][j][2] = dstat[j][2];
			  detmap[tscount][j][3] = dstat[j][3];
			}
		      tscount++;
		      inval++; inval_sub++;
		    }
		}
	      // One pol not consecutive
	      else
		{
		  // Create fake detection with measured mean
		  if(ifverbose)
		    fprintf(stderr,"Gap in frame count detected in file %d subint %d (%f sec). Fake detection with measured mean.\n", pf.filenum, pf.tot_rows, pf.T);
		  getVDIFFrameFakeDetection_1chan(mean,rms,nchan,det,seed,fbytes,dstat);

		  // Copy detection
		  for(j=0;j<nchan;j++)
		    {
		      detmap[tscount][j][0] = dstat[j][0];
		      detmap[tscount][j][1] = dstat[j][1];
		      detmap[tscount][j][2] = dstat[j][2];
		      detmap[tscount][j][3] = dstat[j][3];
		    }
		  tscount++;
		  inval++; inval_sub++;
		}
	    }
	}

      /* ------------- Detection jobs --------------*/
      // Scatter p0, p1 
      MPI_Scatter(bufferall[0], detsum*fbytes, MPI_CHAR, localbuff[0], detsum*fbytes, MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Scatter(bufferall[1], detsum*fbytes, MPI_CHAR, localbuff[1], detsum*fbytes, MPI_CHAR, 0, MPI_COMM_WORLD);
      detarr = (float *)malloc(sizeof(float)*detsum*nchan*4);

      // Get detections
      for(MPIi=0;MPIi<detsum;MPIi++)
	{
	  memcpy(MPIbuff[0], localbuff[0]+MPIi*fbytes, fbytes);
	  memcpy(MPIbuff[1], localbuff[1]+MPIi*fbytes, fbytes);

	  getVDIFFrameDetection_1chan(MPIbuffer[0], MPIbuffer[1], fbytes, MPIdet, nchan, dstat, in_p0, in_p1, out_p0, out_p1, pl0, pl1);

	  // Repack detection into 1D array
	  for(MPIj=0;MPIj<nchan;MPIj++)
	    {
	      detarr[MPIi*nchan*4+MPIj*4+0] = MPIdet[MPIj][0];
	      detarr[MPIi*nchan*4+MPIj*4+1] = MPIdet[MPIj][1];
	      detarr[MPIi*nchan*4+MPIj*4+2] = MPIdet[MPIj][2];
	      detarr[MPIi*nchan*4+MPIj*4+3] = MPIdet[MPIj][3];
	    }
	}

      // Gather detections
      MPI_Gather(detarr, detsum*nchan*4, MPI_FLOAT, detarrglobal, detsum*nchan*4, MPI_FLOAT, 0, MPI_COMM_WORLD);

      /*------------- Writing detections on head node --------------*/
      if(world_rank == 0)
	{

	  // Reassign detection samples
	  for(i=0;i<detsum*nthreads;i++)
	    for(j=0;j<nchan;j++)
	      {
		detmap[valtsmarker[i]][j][0]=detarrglobal[i*nchan*4+j*4+0];
		detmap[valtsmarker[i]][j][0]=detarrglobal[i*nchan*4+j*4+0];
		detmap[valtsmarker[i]][j][0]=detarrglobal[i*nchan*4+j*4+0];
		detmap[valtsmarker[i]][j][0]=detarrglobal[i*nchan*4+j*4+0];
	      }

	  // Loop to offload subints
	  nsuboff = tscount / (pf.hdr.nsblk * tsf);
	  for(k=0; k<nsuboff; k++)
	    {
	      // Initilize subint buffer
	      memset(pf.sub.rawdata, 0, sizeof(unsigned char)*pf.sub.bytes_per_subint);

	      // Fill time samples in subint
	      for(i=0;i<pf.hdr.nsblk;i++)
		{
		  for(j=0;j<nchan;j++)
		    {
		      sdet[j][0]=0.0;
		      sdet[j][1]=0.0;
		      sdet[j][2]=0.0;
		      sdet[j][3]=0.0;
		    }

		  // Create one time sample
		  for(p=0;p<tsf;p++)
		    {
		      for(j=0;j<nchan;j++)
			{
			  sdet[j][0]+=detmap[k*pf.hdr.nsblk*tsf+i*pf.hdr.nsblk+p][j][0];
			  sdet[j][1]+=detmap[k*pf.hdr.nsblk*tsf+i*pf.hdr.nsblk+p][j][1];
			  sdet[j][2]+=detmap[k*pf.hdr.nsblk*tsf+i*pf.hdr.nsblk+p][j][2];
			  sdet[j][3]+=detmap[k*pf.hdr.nsblk*tsf+i*pf.hdr.nsblk+p][j][3];
			}
		    }

		  // Write one time sample
		  for(j=0;j<nchan;j++)
		    {
		      if (npol == 4)
			{
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*j,&sdet[j][0],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*1+sizeof(float)*j,&sdet[j][1],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*2+sizeof(float)*j,&sdet[j][2],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*3+sizeof(float)*j,&sdet[j][3],sizeof(float));
			}
		      else if (npol == 2)
			{
			  memcpy(pf.sub.rawdata+i*sizeof(float)*2*nchan+sizeof(float)*j,&sdet[j][0],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*2*nchan+sizeof(float)*nchan*1+sizeof(float)*j,&sdet[j][1],sizeof(float));
			}
		      else if (npol == 1)
			  memcpy(pf.sub.rawdata+i*sizeof(float)*1*nchan+sizeof(float)*j,&sdet[j][0],sizeof(float));
		    }
		}

	      // Update offset from Start of subint
	      pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;

	      // Write subint
	      psrfits_write_subint(&pf);
	      fprintf(stdout,"Subint written: %d. Faked samples: %lu out of %lu.\n",pf.tot_rows,inval_sub,pf.hdr.nsblk);
	    }

	  // Reloated detection leftovers to the beginning
	  if(tscount % (pf.hdr.nsblk * tsf))
	    {
	      for(i=0;i<tscount-pf.hdr.nsblk*tsf*nsuboff;i++)
		for(j=0;j<nchan;j++)
		  {
		    detmap[i][j][0] = detmap[pf.hdr.nsblk*tsf*nsuboff+i][j][0];
		    detmap[i][j][1] = detmap[pf.hdr.nsblk*tsf*nsuboff+i][j][1];
		    detmap[i][j][2] = detmap[pf.hdr.nsblk*tsf*nsuboff+i][j][2];
		    detmap[i][j][3] = detmap[pf.hdr.nsblk*tsf*nsuboff+i][j][3];
		  }
	      tscount = tscount - pf.hdr.nsblk * tsf * nsuboff;
	    }

	  // Break when subint is not complete
	  if(k!=tsf || i!=pf.hdr.nsblk) break;
	}
    }while(!feof (vdif[0]) && !feof (vdif[1]) && !pf.status && pf.T < pf.hdr.scanlen);
	
  // Close the last file and cleanup

  MPI_Finalize();
  fits_close_file(pf.fptr, &(pf.status));
  free(pf.sub.dat_freqs);
  free(pf.sub.dat_weights);
  free(pf.sub.dat_offsets);
  free(pf.sub.dat_scales);
  free(pf.sub.rawdata);
  free(buffer[0]);
  free(buffer[1]);
  fclose(vdif[0]);
  fclose(vdif[1]);
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);

  printf("Wrote %d subints (%f sec) in %d files.\n",pf.tot_rows, pf.T, pf.filenum);
  printf("Percentage of valid data: %.2f%%\n",(1.0-(float)inval/offset[0])*100.0);

  return;
}
