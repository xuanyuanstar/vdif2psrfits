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
#include "vdifio.h"
#include "vdif2psrfits.h"
#include "psrfits.h"
#include "dec2hms.h"
#include "fftw3.h"

static double VDIF_BW = 2000.0; //Total Bandwidth in MHz
static int VDIF_BIT = 2;   //Bit per sample
static int VDIF_NCHAN = 32; //Number of channels

void usage(char *prg_name)
{
  fprintf(stdout,
		             "%s [options]\n"
		             " -f   Observing central frequency (MHz)\n"
            		 " -i   Input vdif pol0\n"
		             " -j   Input vdif pol1\n"
		             " -n   Band sense (-1 for lower-side band, 1 for upper-side band)\n" 
		             " -s   Seconds to get mean to fill in invalid frames\n"
		             " -t   Time sample scrunch factor (by default 1). One time sample 8 microsecond\n"
		             " -O   Route of the output file(s) \n"
		  " -h   Available options\n",
		  prg_name);
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *vdif[2],*out;

  const vdif_header *header[2];
  struct psrfits pf;
  
  char vname[2][1024], oroute[1024], ut[30],mjd_str[25],vfhdr[2][VDIF_HEADER_BYTES];
  int arg,j_i,j_j,j_O,n_f,mon[2],mon_nxt,i,j,k,t,nfps,fbytes,vd[2],nf_stat,ftot[2][2][VDIF_NCHAN],ct,tsf,bs,tet;
  float freq,s_stat,fmean[2][2][VDIF_NCHAN],dat,stot[2][2][VDIF_NCHAN];
  double spf;
  long double mjd;
  long int idx[2],sec[2],num[2];
  unsigned char *buffer[2];
  int det[VDIF_NCHAN][4];
  float sdet[VDIF_NCHAN][4];
  
  freq=0.0;
  j_i=0;
  j_j=0;
  j_O=0;
  s_stat=0.0;
  tsf=1;
  
  //Read arguments
  while ((arg=getopt(argc,argv,"hf:i:j:s:n:t:O:")) != -1)
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
  
  //Check if arguments are enough to procceed
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
  
  //Read the first header of vdif pol0
  vdif[0]=fopen(vname[0],"rb");
  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
  fclose(vdif[0]);
  header[0]=(const vdif_header *)vfhdr[0];
  
  //Get header info
  /*-----------------------*/
  //Get frame bytes
  fbytes=getVDIFFrameBytes(header[0])-VDIF_HEADER_BYTES;

  //Calculate time interval (in unit of microsecond) of a frame
  spf=(double)fbytes/VDIF_BIT*8/(VDIF_BW*2);
  
  //Calculate MJD
  mjd=(long double)getVDIFFrameMJD(header[0])+(long double)(getVDIFFrameSecond(header[0]))/86400.0;

  //Convert MJD to UT
  mjd2date(mjd,ut);
  printf("UT of start: %s\n",ut);

  //Calculate how many frames to get statistics
  nf_stat=s_stat*1.0e6/spf;
  
  //Allocate memo for frames
  buffer[0]=malloc(sizeof(unsigned char)*fbytes);
  buffer[1]=malloc(sizeof(unsigned char)*fbytes);

  //Scan the beginning specified length of data, choose valid frames to get mean of total value in each frame
  for(j=0;j<2;j++)
	{
	  //Initialize mean value of odd and even samples
	  for(i=0;i<VDIF_NCHAN;i++)
		{
		  fmean[j][0][i]=0.0;
		  fmean[j][1][i]=0.0;
		}
	  ct=0;
	  
	  //Open file
	  vdif[j]=fopen(vname[j],"rb");	
	  for(i=0;i<nf_stat;i++)
		{
		  //Read header
		  fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);
		  header[j]=(const vdif_header *)vfhdr[j];

		  //Valid frame
		  if(!getVDIFFrameInvalid(header[j]))
			{
			  //Read data in frame
			  fread(buffer[j],1,fbytes,vdif[j]);
			  
			  //Get total from odd and even samples
			  getVDIFFrameTotal_32chan(buffer[j],fbytes,&ftot[j]);

			  //Add to the mean calulation
			  for(k=0;k<VDIF_NCHAN;k++)
				{
				  fmean[j][0][k]+=(float)ftot[j][0][k];
				  fmean[j][1][k]+=(float)ftot[j][1][k];
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
	  for(k=0;k<VDIF_NCHAN;k++)
		{
		  fmean[j][0][k]/=ct;
		  fmean[j][1][k]/=ct;
		  printf("Mean in chan%i pol%i: %f  %f\n",k,j,fmean[j][0][k],fmean[j][1][k]);
		}	  
	}

  // Open VDIF files for reading
  for(j=0;j<2;j++)
	vdif[j]=fopen(vname[j],"rb");
  
  // Calibrate the difference in starting time between the two pols
  printf("Calibrating potential difference in starting time between two pols...\n");
  for(j=0;j<2;j++)
	{
	  fread(vfhdr[j],1,VDIF_HEADER_BYTES,vdif[j]);
	  header[j]=(const vdif_header *)vfhdr[j];
	  
	  //Get epoch in unit of 6-month period
	  mon[j]=getVDIFFrameMJD(header[j]);

	  //Get second after epoch
	  sec[j]=getVDIFFrameSecond(header[j]);

	  //Get number of frame after second
	  num[j]=getVDIFFrameNumber(header[j]);
	}  
  while(mon[0]!=mon[1])
	{
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

  mjd=(long double)mon[0]+(long double)sec[0]/86400.0+(long double)spf*num[0]/86400.0;

  //Set psrfits main header
  printf("Setting up PSRFITS output...\n");
  pf.filenum = 0;           // This is the crucial one to set to initialize things
  pf.rows_per_file = 200;  // Need to set this based on PSRFITS_MAXFILELEN

  //Set values for our hdrinfo structure
  pf.hdr.scanlen = 86400; // in sec
  strcpy(pf.hdr.observer, "A. Eintein");
  strcpy(pf.hdr.telescope, "ALMA");
  strcpy(pf.hdr.obs_mode, "SEARCH");
  strcpy(pf.hdr.backend, "MARK6");
  strcpy(pf.hdr.source, "J0835-4510");
  strcpy(pf.hdr.date_obs, ut);
  strcpy(pf.hdr.poln_type, "LIN");
  strcpy(pf.hdr.poln_order, "AABBCRCI");
  strcpy(pf.hdr.track_mode, "TRACK");
  strcpy(pf.hdr.cal_mode, "OFF");
  strcpy(pf.hdr.feed_mode, "FA");
  pf.hdr.dt = spf*tsf/1.0e6;
  pf.hdr.fctr = freq;
  pf.hdr.BW = VDIF_BW*bs;
  pf.hdr.nchan = VDIF_NCHAN;
  pf.hdr.MJD_epoch = mjd;
  pf.hdr.ra2000 = 302.0876876;
  dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
  pf.hdr.dec2000 = -3.456987698;
  dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);
  pf.hdr.azimuth = 123.123;
  pf.hdr.zenith_ang = 23.0;
  pf.hdr.beam_FWHM = 0.25;
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
  pf.hdr.npol = 4;
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
  
  //Set values for our subint structure
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

  printf("Header prepared. Start to write data...\n");
  
  // Main loop to write subints
  do
	{
	  // Fill time samples in each subint: pf.sub.rawdata
	  for(i=0;i<pf.hdr.nsblk;i++)
		{
		  for(k=0;k<VDIF_NCHAN;k++)
			{
			  sdet[k][0]=0.0;
			  sdet[k][1]=0.0;
			  sdet[k][2]=0.0;
			  sdet[k][3]=0.0;
			} 
		  // Loop over frames
		  for(k=0;k<tsf;k++)
			{
			  //Read frame header
			  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
			  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
			  header[0]=(const vdif_header *)vfhdr[0];
			  header[1]=(const vdif_header *)vfhdr[1];

			  //Valid frame
			  if(!getVDIFFrameInvalid(header[0]) && !getVDIFFrameInvalid(header[1]))
				{
				  //Read data in frame
				  fread(buffer[0],1,fbytes,vdif[0]);
				  fread(buffer[1],1,fbytes,vdif[1]);

				  getVDIFFrameDetection_coherence_32chan(buffer[0],buffer[1],fbytes,det);
				  
				  //Accumulate value
				  for(j=0;j<VDIF_NCHAN;j++)
					{
					  sdet[j][0]+=(float)det[j][0];
					  sdet[j][1]+=(float)det[j][1];
					  sdet[j][2]+=(float)det[j][2];
					  sdet[j][3]+=(float)det[j][3];
					}
				}
			  //Invalide frame
			  else
				{
				  //Skip the data
				  fseek(vdif[0],fbytes,SEEK_CUR);
				  fseek(vdif[1],fbytes,SEEK_CUR);

				  //Use mean of total
				  for(j=0;j<VDIF_NCHAN;j++)
					{
					  sdet[j][0]+=(float)det[j][0];
					  sdet[j][1]+=(float)det[j][1];
					  sdet[j][2]+=(float)det[j][2];
					  sdet[j][3]+=(float)det[j][3];
					}
				  fprintf(stderr,"Invalid frame detected. Supply both pols with mean of total.\n",j);
				}			  
			  //Break at the end of vdif file
			  if(feof (vdif[0]) || feof (vdif[1])) break;
			}		  
		  //Break when not enough frames to get a sample
		  if(k!=tsf) break;
		  
		  //Write detections in pf.sub.rawdata, in 32-bit float and FPT order (freq, pol, time)
		  for(j=0;j<VDIF_NCHAN;j++)
			{
			  //printf("%f %f %f %f\n",sdet[j][0],sdet[j][1],sdet[j][2],sdet[j][3]);
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN+sizeof(float)*j,&sdet[j][1],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*2+sizeof(float)*j,&sdet[j][2],sizeof(float));
			  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*3+sizeof(float)*j,&sdet[j][3],sizeof(float));
			}
		}

	  //Break when subint is not complete
	  if(k!=tsf || i!=pf.hdr.nsblk) break;

	  //Update offset from Start of subint
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;

	  //Write subint
	  psrfits_write_subint(&pf);
	  printf("Subint %i written.\n",pf.sub.tsubint);
	  
	} while(!feof (vdif[0]) && !feof (vdif[1]) && !pf.status && pf.T < pf.hdr.scanlen);
	
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

  return 0;
}
