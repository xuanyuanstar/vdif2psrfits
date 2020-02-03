//Convert vdif format to psrfits search mode format in 32-bit float
//Input vdif configuration is 2-bit real-sampled, 2048 MHz bandwdith, one frequency channel, two pols in separated files
//Sampling interval 1/2048/2 microsecond
//One frame 8192 bytes, thus 8 microsecond per frame (one frame contains only one pol)


#define VDIF_HEADER_BYTES       32

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "vdifio.h"
#include "psrfits.h"
#include "vdif2psrfits.h"
#include "dec2hms.h"
#include "ran.c"

static double VDIF_BW = 2000.0; //Total Bandwidth in MHz
static int VDIF_BIT = 2;   //Bit per sample
static int VDIF_NCHAN = 32; //Number of channels

void usage(char *prg_name)
{
  fprintf(stdout,
		  "%s [options]\n"
		  "  -f      Observing central frequency (MHz)\n"
                  "  -i      Input vdif pol0\n"
		  "  -j      Input vdif pol1\n"
		  "  -n      Band sense (-1 for lower-side, 1 for upper-side, by default 1)\n" 
		  "  -s      Seconds to get statistics to fill in invalid frames\n"
                  "  -k      Seconds to skip from the beginning (default 0)\n"
		  "  -t      Time sample scrunch factor (by default 1). One time sample 8 microsecond\n"
                  "  -S      Name of the source (by default J0000+0000)\n"
	          "  -r      RA of the source (AA:BB:CC.DD)\n"
                  "  -c      Dec of the source (+AA:BB:CC.DD)\n"
		  "  -D      Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, S for Stokes, P for polarised signal, S for stokes, by default C)\n"
		  "  -O      Route of the output file(s).\n"
		  "  -h      Available options\n"
		  "\n"
		  "Patching power dip options:\n"
		  "  -P               Replace power dip raw samples with random noise \n"
		  "  -M               Replace power dip detections with mean \n"
		  "  -p               Starting phase of data in scan+dip cycle (by default 0)\n",
		  prg_name);
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *vdif[2],*out;

  const vdif_header *header[2];
  struct psrfits pf;
  
  char vname[2][1024], oroute[1024], ut[30],mjd_str[25],vfhdr[2][VDIF_HEADER_BYTES],srcname[16],dstat,ra[64],dec[64];
  int arg,j_i,j_j,j_O,n_f,mon[2],mon_nxt,i,j,k,p,nfps,fbytes,fnum,vd[2],nf_stat,ftot[2][2][VDIF_NCHAN],ct,tsf,bs,tet,nf_skip,dati,npol,pch,mean_sampl,sk,ftotct;
  float freq,s_stat,dat,s_skip;
  double spf,pha_start,len_scan,len_dip,mean_det[VDIF_NCHAN][4],acc_det[VDIF_NCHAN][4],rms_det[VDIF_NCHAN][4],accsq_det[VDIF_NCHAN][4];
  long double mjd;
  long int idx[2],sec[2],num[2],seed,iseed,pha_start_nf,len_scan_nf,len_dip_nf,pha_ct,fct;
  unsigned char *buffer[2], *obuffer[2];
  float det[VDIF_NCHAN][4],sdet[VDIF_NCHAN][4];
  time_t t;
  
  freq=0.0;
  j_i=0;
  j_j=0;
  j_O=0;
  s_stat=0.0;
  tsf=1;
  s_skip=0.0;
  strcpy(srcname,"J0000+0000");
  dstat='C';
  npol=4;
  pch=0;
  pha_start=0.0;
  len_scan=16.128;
  len_dip=2.064;
  mean_sampl=0;
  ftotct=0;

  if(argc==1)
    {
      usage(argv[0]);
      exit(0);
    }
  
  // Read arguments
  while ((arg=getopt(argc,argv,"hf:i:j:s:n:k:t:O:S:D:r:c:dPp:M")) != -1)
	{
	  switch(arg)
		{
		case 'f':
		  freq=atof(optarg);
		  break;
		  
		case 'i':
		  strcpy(vname[0],optarg);
		  j_i=1;
		  break;

		case 'j':
		  strcpy(vname[1],optarg);
		  j_j=1;
		  break;

		case 's':
		  s_stat=atof(optarg);
		  break;

		case 'n':
		  bs=atoi(optarg);
		  break;
		  
		case 't':
		  tsf=atoi(optarg);
		  break;

		case 'k':
		  s_skip=atof(optarg);
		  break;

		case 'S':
		  strcpy(srcname,optarg);
		  break;

		case 'r':
		  strcpy(ra,optarg);
		  break;

		case 'c':
		  strcpy(dec,optarg);
		  break;
		  
		case 'D':
		  strcpy(&dstat,optarg);
		  if(dstat=='C' || dstat=='S')
			npol=4;
		  else if(dstat=='P')
			npol=2;
		  else
			npol=1;
		  break;
		  
		case 'P':
		  pch=1;
		  break;

		case 'M':
		  pch=2;
		  break;
		  
		case 'p':
		  pha_start=atof(optarg);
		  break;
		  
		case 'O':
		  strcpy(oroute,optarg);
		  j_O=1;
		  break;

		case 'h':
		  usage(argv[0]);
		  return 0;

		default:
		  usage(argv[0]);
		  return 0;
		}
	}
  
  // Check if arguments are enough to procceed
  if(freq==0.0)
	{
	  fprintf(stderr,"Missing info of observing central frequency.\n");
	  exit(0);
	}
  if(j_i==0)
	{
	  fprintf(stderr,"No input file provided for pol0.\n");
	  exit(0);
	}
  if(j_j==0)
	{
	  fprintf(stderr,"No input file provided for pol1.\n");
	  exit(0);
	}
  if(tsf<1)
	{
	  fprintf(stderr,"Invalid sample scrunching factor.\n");
	  exit(0);
	}  
  if(j_O==0)
	{
	  fprintf(stderr,"No output route specified.\n");
	  exit(0);
	}
  if(bs!=-1 && bs!=1)
	{
	  fprintf(stderr,"Not readable band sense.\n");
	  exit(0);
	}
  if(dstat!='I' && dstat!='C' && dstat!='X' && dstat!='Y' && dstat!='S' && dstat!='P')
	{
	  fprintf(stderr,"Not recognized status for output data.\n");
	  exit(0);
	}
  
  // Get seed for random generator
  time(&t);
  iseed=0-t;
  seed=iseed;
  
  // Read the first header of vdif pol0
  vdif[0]=fopen(vname[0],"rb");
  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
  fclose(vdif[0]);
  header[0]=(const vdif_header *)vfhdr[0];
  
  // Get header info
  /*-----------------------*/
  // Get frame bytes
  fbytes=getVDIFFrameBytes(header[0])-VDIF_HEADER_BYTES;

  // Calculate time interval (in unit of microsecond) of a frame
  spf=(double)fbytes/VDIF_BIT*8/(VDIF_BW*2);
  
  // Calculate how many frames to get statistics
  nf_stat=s_stat*1.0e6/spf;

  // Calculate how many frames to skip from the beginning
  nf_skip=s_skip*1.0e6/spf;
  printf("Number of frames to skip from the beginning: %i.\n",nf_skip);
    
  // Allocate memo for frames
  for(j=0;j<2;j++)
	{
	  buffer[j]=malloc(sizeof(unsigned char)*fbytes);
	  obuffer[j]=malloc(sizeof(unsigned char)*fbytes*4);
	}
  
  // Scan the beginning specified length of data, choose valid frames to get mean of total value in each frame
  fprintf(stderr,"Scan %.2f s data to get statistics, after skipping the first %.2f data...\n",s_stat,s_skip);
  for(j=0;j<4;j++)
	for(k=0;k<VDIF_NCHAN;k++)
	  {
		acc_det[k][j]=0.0;
		accsq_det[k][j]=0.0;
	  }
  fct=0;	  
  for(j=0;j<2;j++)
	{
	  vdif[j]=fopen(vname[j],"rb");

	  // Skip the first given length of data
	  for(i=0;i<nf_skip;i++)
		fseek(vdif[j],VDIF_HEADER_BYTES+fbytes,SEEK_CUR);
	}
  for(i=0;i<nf_stat;i++)
	{
	  // Read header
	  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
	  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
	  header[0]=(const vdif_header *)vfhdr[0];
	  header[1]=(const vdif_header *)vfhdr[1];
		  
	  // Valid frame
	  if(!getVDIFFrameInvalid(header[0]) && !getVDIFFrameInvalid(header[1]))
		{
		  // Read data in frame
		  fread(buffer[0],1,fbytes,vdif[0]);
		  fread(buffer[1],1,fbytes,vdif[1]);
		  
		  // Accumulate values for detection mean
		  getVDIFFrameDetection_32chan(buffer[0],buffer[1],fbytes,det,dstat);
		  for(j=0;j<VDIF_NCHAN;j++)
			for(p=0;p<4;p++)
			  {
				acc_det[j][p]+=(double)det[j][p];
				accsq_det[j][p]+=pow((double)det[j][p],2.0);
			  }
		  fct++;
		}
	  // Invalid frame
	  else
		{
		  fseek(vdif[0],fbytes,SEEK_CUR);
		  fseek(vdif[1],fbytes,SEEK_CUR);
		}
	}
  fclose(vdif[0]);
  fclose(vdif[1]);

  // Initialize patching param.
  for(j=0;j<VDIF_NCHAN;j++)
	for(p=0;p<4;p++)
	  {
		mean_det[j][p]=acc_det[j][p]/(double)fct;
		rms_det[j][p]=sqrt(accsq_det[j][p]/fct-pow(mean_det[j][p],2.0));
		fprintf(stderr,"Mean & rms det chan%i, pol%i: %lf %lf\n",j,p,mean_det[j][p],rms_det[j][p]);
	  }
  
  // Open VDIF files and skip the first given length of data
  for(j=0;j<2;j++)
	{
	  vdif[j]=fopen(vname[j],"rb");
	  for(i=0;i<nf_skip;i++)
		fseek(vdif[j],VDIF_HEADER_BYTES+fbytes,SEEK_CUR);
	}
  
  // Calibrate the difference in starting time between the two pols
  printf("Calibrating potential difference in starting time between two pols...\n");
  for(j=0;j<2;j++)
	{
	  fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);
	  header[j]=(const vdif_header *)vfhdr[j];
	  
	  // Get epoch in unit of 6-month period
	  mon[j]=getVDIFFrameMJD(header[j]);

	  // Get second after epoch
	  sec[j]=getVDIFFrameSecond(header[j]);

	  // Get number of frame after second
	  num[j]=getVDIFFrameNumber(header[j]);
	}  
  while(mon[0]!=mon[1])
	{
	  printf("Mon not same %i %i\n",mon[0],mon[1]);
	  if(mon[0]<mon[1])
		{
		  fseek(vdif[0],fbytes,SEEK_CUR);
		  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
		  header[0]=(const vdif_header *)vfhdr[0];
		  mon[0]=getVDIFFrameMJD(header[0]);
		  sec[0]=getVDIFFrameSecond(header[0]);
		}
	  else
		{
		  fseek(vdif[1],fbytes,SEEK_CUR);
		  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
		  header[1]=(const vdif_header *)vfhdr[1];
		  mon[1]=getVDIFFrameMJD(header[1]);
		  sec[1]=getVDIFFrameSecond(header[1]);
		}
	}
  while(sec[0]!=sec[1])
	{
	  printf("sec not same\n");
	  if(sec[0]<sec[1])
		{
		  fseek(vdif[0],fbytes,SEEK_CUR);
		  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
		  header[0]=(const vdif_header *)vfhdr[0];
		  sec[0]=getVDIFFrameSecond(header[0]);
		  num[0]=getVDIFFrameNumber(header[0]);
		}
	  else
		{
		  fseek(vdif[1],fbytes,SEEK_CUR);
		  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
		  header[1]=(const vdif_header *)vfhdr[1];
		  sec[1]=getVDIFFrameSecond(header[1]);
		  num[1]=getVDIFFrameNumber(header[1]);
		}
	}
  while(num[0]!=num[1])
	{
	  printf("num not same\n");
	  if(num[0]<num[1])
		{
		  fseek(vdif[0],fbytes,SEEK_CUR);
		  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
		  header[0]=(const vdif_header *)vfhdr[0];
		  num[0]=getVDIFFrameNumber(header[0]);
		}
	  else
		{
		  fseek(vdif[1],fbytes,SEEK_CUR);
		  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
		  header[1]=(const vdif_header *)vfhdr[1];
		  num[1]=getVDIFFrameNumber(header[1]);
		}
	}
  for(j=0;j<2;j++)
	fseek(vdif[j],-VDIF_HEADER_BYTES,SEEK_CUR);
  printf("Starting time synchronized. MJD: %i; Second: %i; Number: %i.\n",mon[0],sec[0],num[0]);

  // Get starting MJD and UT
  mjd=(long double)mon[0]+(long double)sec[0]/86400.0+(long double)spf*num[0]/86400.0*1.0e-6;
  mjd2date(mjd,ut);
  printf("start UT of VDIF: %s\n",ut);
  
  // Set psrfits main header
  printf("Setting up PSRFITS output...\n");
  pf.filenum = 0;           // This is the crucial one to set to initialize things
  pf.rows_per_file = 200;  // Need to set this based on PSRFITS_MAXFILELEN

  // Set values for our hdrinfo structure
  pf.hdr.scanlen = 86400; // in sec
  strcpy(pf.hdr.observer, "A. Eintein");
  strcpy(pf.hdr.telescope, "ALMA");
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
  pf.hdr.fctr = freq;
  pf.hdr.BW = VDIF_BW*bs;
  pf.hdr.nchan = VDIF_NCHAN;
  pf.hdr.MJD_epoch = mjd;
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
  pf.hdr.nsblk = 8192;
  pf.hdr.ds_time_fact = 1;
  pf.hdr.ds_freq_fact = 1;
  sprintf(pf.basefilename, "%s/%s",oroute,ut);
  
  psrfits_create(&pf);
  
  // Set values for our subint structure
  pf.sub.tsubint = pf.hdr.nsblk * pf.hdr.dt;
  pf.tot_rows = 0.0;
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

  // Initialize param. for patching
  pha_start_nf = lround(pha_start * (len_scan + len_dip) / (pf.hdr.dt/tsf) );
  len_scan_nf = len_scan / (pf.hdr.dt/tsf);
  len_dip_nf = len_dip / (pf.hdr.dt/tsf);
  pha_ct = pha_start_nf;
  fct=0;
  for(j=0;j<4;j++)
	for(k=0;k<VDIF_NCHAN;k++)
	  {
		acc_det[k][j]=0.0;
		accsq_det[k][j]=0.0;
	  }
  
  fprintf(stderr,"Header prepared. Start to write data...\n");

  // Main loop to write subints
  do
	{
	  // Fill time samples in each subint: pf.sub.rawdata
	  for(i=0;i<pf.hdr.nsblk;i++)
		{
		  // Initialize detection block
		  for(k=0;k<VDIF_NCHAN;k++)
			{
			  for(p=0;p<npol;p++)
				sdet[k][p]=0.0;
			} 
		  // Loop over frames
		  for(k=0;k<tsf;k++)
			{
			  // Read frame header
			  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
			  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
			  header[0]=(const vdif_header *)vfhdr[0];
			  header[1]=(const vdif_header *)vfhdr[1];
			  ftotct++;
			  
			  // Check if a proper frame
			  fbytes=getVDIFFrameBytes(header[0])-VDIF_HEADER_BYTES;
			  num[0]=getVDIFFrameNumber(header[0]);
			  num[1]=getVDIFFrameNumber(header[1]);
			  if(fbytes < 32 || fbytes > 144512) {
				printf("False number of framesize in pol0 at frame %i.\n",ftotct);
				exit(0);
			  }
			  fbytes=getVDIFFrameBytes(header[1])-VDIF_HEADER_BYTES;
			  if(fbytes < 32 || fbytes > 144512) {
				printf("False number of framesize in pol1 at frame %i.\n",ftotct);
				exit(0);
			  }
			  if(num[0] != num[1]) printf("Warning: frame number start to deviate.\n");
			  
			  // Valid frame
			  if(!getVDIFFrameInvalid(header[0]) && !getVDIFFrameInvalid(header[1]))
				{
				  // Read data in frame
				  fread(buffer[0],1,fbytes,vdif[0]);
				  fread(buffer[1],1,fbytes,vdif[1]);

				  // Get detection
				  getVDIFFrameDetection_32chan(buffer[0],buffer[1],fbytes,det,dstat);

				  // Subscan phase
				  if(pha_ct < len_scan_nf)
					{				  
					  // Accumulate values for detection mean
					  for(j=0;j<VDIF_NCHAN;j++)
						for(p=0;p<4;p++)
						  {
							acc_det[j][p]+=(double)det[j][p];
							accsq_det[j][p]+=pow((double)det[j][p],2.0);
						  }
					  fct++;
					}
				  else // Dip phase
					{
					  // Update patching param, when entering dip phase
					  if(pha_ct == len_scan_nf)
						{
						  // Update detection mean
						  for(j=0;j<VDIF_NCHAN;j++)
							for(p=0;p<4;p++)
							  {
								// In case no values in the past cycle
								if(acc_det[j][p]!=0.0)
								  {
									mean_det[j][p]=acc_det[j][p]/(double)fct;
									rms_det[j][p]=sqrt(accsq_det[j][p]/(double)fct-pow(mean_det[j][p],2.0));
								  }
								acc_det[j][p]=0.0;
								accsq_det[j][p]=0.0;
							  }
						  fct=0;
						}
                      // Patch detection with mean+rms
					  if(pch == 1)
						{
						  for(j=0;j<VDIF_NCHAN;j++)
							for(p=0;p<4;p++)
							  det[j][p]=mean_det[j][p]+rms_det[j][p]*gasdev(&seed);
						}
					  else if(pch == 2)
						{
						  // Patch fake detection from mean
						  for(j=0;j<VDIF_NCHAN;j++)
							for(p=0;p<4;p++)
							  det[j][p]=mean_det[j][p];
						}
					}
				}
			  else // Invalid frame	  
				{
				  // Fake detection from mean
				  fprintf(stderr,"Invalid frame detected. Use measured mean to generate fake detection.\n");
				  for(j=0;j<VDIF_NCHAN;j++)
					for(p=0;p<4;p++)
					  det[j][p]=mean_det[j][p];

				  // Skip the data
				  fseek(vdif[0],fbytes,SEEK_CUR);
				  fseek(vdif[1],fbytes,SEEK_CUR);
				}
			
			  // Accumulate detection value
			  for(j=0;j<VDIF_NCHAN;j++)
				{
				  for(p=0;p<npol;p++)
					sdet[j][p]+=det[j][p];
				}

			  // Update phase counters
			  pha_ct++;
				  
			  // Reset phase and counter
			  if(pha_ct == (len_scan_nf + len_dip_nf))
				{
				  pha_ct = 0;

				  // Start another sequence for ran
				  iseed=iseed-2;
				  seed=iseed;
				}

			  // Break at the end of vdif file
			  if(feof (vdif[0]) || feof (vdif[1])) break;
			}
		
		  // Break when not enough frames to get a sample
		  if(k!=tsf) break;
		  
		  // Write detections in pf.sub.rawdata, in 32-bit float and FPT order (freq, pol, time)
		  for(j=0;j<VDIF_NCHAN;j++)
			{
			  if (npol == 4)
				{
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*1+sizeof(float)*j,&sdet[j][1],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*2+sizeof(float)*j,&sdet[j][2],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*3+sizeof(float)*j,&sdet[j][3],sizeof(float));
				}
			  else if (npol == 2)
				{
                  memcpy(pf.sub.rawdata+i*sizeof(float)*2*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*2*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*1+sizeof(float)*j,&sdet[j][1],sizeof(float));
				}
			  else if (npol == 1)
				{
				  memcpy(pf.sub.rawdata+i*sizeof(float)*1*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
				}
			}
		  // Break at the end of vdif file
		  if(feof (vdif[0]) || feof (vdif[1])) break;
		}

	  // Break when subint is not complete
	  if(k!=tsf || i!=pf.hdr.nsblk) break;

	  // Update offset from Start of subint
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;

	  // Write subint
	  psrfits_write_subint(&pf);
	  printf("Subint %i written.\n",pf.sub.tsubint);
	  if(feof(vdif[0]) || feof(vdif[1])) break;
	} while(!feof(vdif[0]) && !feof(vdif[1]) && !pf.status && pf.T < pf.hdr.scanlen);
	
  // Close the last file and cleanup
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
  
  printf("Wrote %d subints (%f sec) in %d files.\n",pf.tot_rows, pf.T, pf.filenum);

  return;
}
