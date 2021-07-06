#include "cvrt2to8.c"
#include "ran.c"
#include <malloc.h>
#include <complex.h>
#include <stdbool.h>
#include <fftw3.h>
#include "vdifio.h"

// Mean of unsigned 2-bit samples
static float mean2bspl = 1.5;

void getDetection(float p0r, float p0i, float p1r, float p1i, float *det, char dstat)
{
  if (dstat == 'C')
	{
	  det[0]=p0r*p0r+p0i*p0i;
	  det[1]=p1r*p1r+p1i*p1i;
	  det[2]=p0r*p1r+p0i*p1i;
	  det[3]=p0r*p1i-p0i*p1r;
	}
  else if (dstat == 'I')
	det[0]=p0r*p0r+p0i*p0i+p1r*p1r+p1i*p1i;
  else if (dstat == 'X')
	det[0]=p0r*p0r+p0i*p0i;
  else if (dstat == 'Y')
	det[0]=p1r*p1r+p1i*p1i;
  else if (dstat == 'S')
	{
	  det[0]=(p0r*p0r+p0i*p0i)+(p1r*p1r+p1i*p1i);
	  det[1]=(p0r*p0r+p0i*p0i)-(p1r*p1r+p1i*p1i);
	  det[2]=2.0*(p0r*p1r+p0i*p1i);
	  det[3]=2.0*(p0r*p1i-p0i*p1r); // PSR/IEEE convention in Stokes V
	}
  else if (dstat == 'P')
	{
	  det[0]=sqrt(pow((p0r*p0r+p0i*p0i)-(p1r*p1r+p1i*p1i),2.0)+pow(2.0*(p0r*p1r+p0i*p1i),2.0));
	  det[1]=2.0*(p0r*p1i-p0i*p1r);
	}
}

// Get coherence detection from real 2-bit, 32-channel frame of two pols
void getVDIFFrameDetection_32chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, float det[][4],char dstat, float *in_p0, float *in_p1, fftwf_complex *out_p0, fftwf_complex *out_p1, fftwf_plan pl0, fftwf_plan pl1)
{
  float dets[4];
  int i,j,k,Nts,Nchan,npol;
  unsigned char *buff8[2];

  // Decode dstat to get npol
  if(dstat=='C' || dstat=='S')
	npol=4;
  else if(dstat=='P')
	npol=2;
  else
	npol=1;

  Nchan=32;
  for(i=0;i<2;i++)
    buff8[i]=malloc(sizeof(unsigned char)*fbytes*4);

  //Extend from 2-bit to 8-bit
  convert2to8(buff8[0], src_p0, fbytes);
  //for(i=0;i<fbytes;i++)
  //{
  //  buff8[0][4*fbytes] = (src_p0[i] & 192) >> 0x6 ;
  //  buff8[1][4*fbytes] = (src_p1[i] & 192) >> 0x6 ;
  ///  buff8[0][4*fbytes+1] = (src_p0[i] & 48) >> 0x4 ;
  // buff8[1][4*fbytes+1] = (src_p1[i] & 48) >> 0x4 ;
  //  buff8[0][4*fbytes+2] = (src_p0[i] & 12) >> 0x2 ;
  //  buff8[1][4*fbytes+2] = (src_p1[i] & 12) >> 0x2 ;
  //  buff8[0][4*fbytes+3] = (src_p0[i] & 3) >> 0x0 ;
  //  buff8[1][4*fbytes+3] = (src_p1[i] & 3) >> 0x0 ;
  //}
  convert2to8(buff8[1], src_p1, fbytes);

  //Number of time samples
  Nts=fbytes*4/Nchan;
  
  //Initialize
  for(k=0;k<Nchan;k++)
	for(j=0;j<4;j++)
	  det[k][j]=0.0;
  for(j=0;j<4;j++)
    dets[j]=0.0;
  
  //Detect each channel
  for(k=0;k<Nchan;k++)
    {
      //Read time series
      for(i=0;i<Nts;i++)
	{
	  if(k<Nchan/2)
	    {
	      //in_p0[i]=1.0;in_p1[i]=1.0;
	      in_p0[i]=(float)((int)buff8[0][i*Nchan+15-k])-mean2bspl;
	      in_p1[i]=(float)((int)buff8[1][i*Nchan+15-k])-mean2bspl;
	    }
	  else
	    {
	      //in_p0[i]=1.0;in_p1[i]=1.0;
	      in_p0[i]=(float)((int)buff8[0][i*Nchan+15+Nchan-k])-mean2bspl;
	      in_p1[i]=(float)((int)buff8[1][i*Nchan+15+Nchan-k])-mean2bspl;
	    }
	}

      // Launch FFT
      fftwf_execute(pl0);
      fftwf_execute(pl1);

      // Make detection for each FFT channel and sum up
      for(i=0;i<=Nts/2;i++)
	{
	  getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);
	  for(j=0;j<npol;j++)
	    det[k][j]+=dets[j];
	}
    }

  //Free up memo
  for(i=0;i<2;i++)
    free(buff8[i]);
}


// Generate fake detection with given mean and rms for a frame with 32 channels
void getVDIFFrameFakeDetection_32chan(double mean_scan[][32], double rms_scan[][32], float det[][4], long int *seed, int fbytes, char dstat)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4], rum;
  int i,j,k,Nts,npol,idat,Nchan;

  Nchan=32;
  
  // Decode dstat to get npol
  if(dstat=='C' || dstat=='S')
	npol=4;
  else if(dstat=='P')
	npol=2;
  else
	npol=1;

  // Number of time samples
  Nts=fbytes*4/Nchan;
  
  // Memo for FFT for two pols
  in = (float **)malloc(sizeof(float *)*2);
  for(i=0;i<2;i++)
	in[i]= (float *)malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));

  // Initialize
  for(k=0;k<Nchan;k++)
	for(j=0;j<4;j++)
	  det[k][j]=0.0;

  for(j=0;j<4;j++)
	dets[j]=0.0;

  // Detect each channel
  for(k=0;k<Nchan;k++)
	{
	  // Operate each pol
	  for(j=0;j<2;j++)
		{
		  // Generate time series with given mean & rms
		  for(i=0;i<Nts;i++)
			in[j][i]=mean_scan[j][k]-mean2bspl;
		}

      // Make FFT plan and execute
	  if(k==0)
		{
		  pl0 = fftwf_plan_dft_r2c_1d(Nts, in[0], out_p0, FFTW_ESTIMATE);
		  pl1 = fftwf_plan_dft_r2c_1d(Nts, in[1], out_p1, FFTW_ESTIMATE);
		}
	  fftwf_execute(pl0);
	  fftwf_execute(pl1);

      // Make detection for each FFT channel and sum up
	  for(i=1;i<Nts/2;i++)
		{
		  getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);
		  for(j=0;j<npol;j++)
			det[k][j]+=dets[j];
		}
	  // Add DC and Nyquist power
	  getDetection(creal(out_p0[0]),cimag(out_p0[0]),creal(out_p1[0]),cimag(out_p1[0]),dets,dstat);
	  for(j=0;j<npol;j++)
		det[k][j]+=dets[j]/2;
	  getDetection(creal(out_p0[Nts/2]),cimag(out_p0[Nts/2]),creal(out_p1[Nts/2]),cimag(out_p1[Nts/2]),dets,dstat);
	  for(j=0;j<npol;j++)
		det[k][j]+=dets[j]/2;
	}

  //Free up memo
  for(i=0;i<2;i++)
	free(in[i]);
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);
}

// Generate fake detection with given mean and rms for a frame with one channel
void getVDIFFrameFakeDetection_1chan(double *mean, double *rms, int Nchan, float det[][4], long int *seed, int fbytes, char dstat)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts,npol, chw;

  // Decode dstat to get npol
  if(dstat=='C' || dstat=='S')
	npol=4;
  else if(dstat=='P')
	npol=2;
  else
	npol=1;
  
  // Number of time samples
  Nts = fbytes*4;

  //Number of FFT spectral per channel
  chw = Nts/2/Nchan;
  
  // Memo for FFT for two pols
  in = (float **)malloc(sizeof(float *)*2);
  for(i=0;i<2;i++)
    in[i]= (float *)malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));

  // Initialize
  for(k=0;k<Nchan;k++)
    for(j=0;j<4;j++)
      det[k][j]=0.0;

  for(j=0;j<4;j++)
    dets[j]=0.0;

  // Operate each pol
  for(j=0;j<2;j++)
    {
      // Generate time series with given mean & rms
      for(i=0;i<Nts;i++)
	//in[j][i]=mean[j]+gasdev(seed)*rms[j]-mean2bspl;
	in[j][i]=mean[j]-mean2bspl;
    }

  // Perform FFT
  pl0 = fftwf_plan_dft_r2c_1d(Nts, in[0], out_p0, FFTW_ESTIMATE);
  pl1 = fftwf_plan_dft_r2c_1d(Nts, in[1], out_p1, FFTW_ESTIMATE);
  fftwf_execute(pl0);
  fftwf_execute(pl1);

  // Make detection
  for(i=1;i<=Nts/2;i++)
    {
      getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);

      // Channel index
      j=(i-1)/chw;

      // Envalue
      for(k=0;k<4;k++)
	det[j][k]+=dets[k];
    }

  //Free up memo
  for(i=0;i<2;i++)
    free(in[i]);
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);
}

void getVDIFFrameDetection_1chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, float det[][4], int nchan, char dstat, float *in_p0, float *in_p1, fftwf_complex *out_p0, fftwf_complex *out_p1, fftwf_plan pl0, fftwf_plan pl1)
{
  float dets[4];
  int i,j,k,Nts,chw;
  unsigned char *buff8[2];
  
  for(i=0;i<2;i++)
    buff8[i]=malloc(sizeof(unsigned char)*fbytes*4);

  //Extend from 2-bit to 8-bit
  convert2to8(buff8[0], src_p0, fbytes);
  convert2to8(buff8[1], src_p1, fbytes);

  //Number of time samples
  Nts=fbytes*4;

  //Number of FFT spectral per channel
  chw=Nts/2/nchan;
  
  // Initialization
  for(j=0;j<4;j++)
    dets[j]=0.0;
  for(j=0;j<nchan;j++)
    for(k=0;k<4;k++)
      det[j][k]=0.0;

  for(i=0;i<Nts;i++)
    {
      in_p0[i] = (float)((int)buff8[0][i])-mean2bspl;
      in_p1[i] = (float)((int)buff8[1][i])-mean2bspl;
    }

  fftwf_execute(pl0);
  fftwf_execute(pl1);

  //Make detection for each FFT channel and sum up to given nchan
  for(i=1;i<=Nts/2;i++)
    {
      getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);

      // Channel index
      j=(i-1)/chw;

      // Envalue
      for(k=0;k<4;k++)
	det[j][k]+=dets[k];
    }

  //Free up memo
  for(i=0;i<2;i++)
    free(buff8[i]);
}

int getVDIFFrameInvalid_robust(const vdif_header *header, int framebytes, bool ifverbose)
{
  int f, mjd,sec,num,inval;
  f=getVDIFFrameBytes(header);
  mjd=getVDIFFrameMJD(header);
  sec=getVDIFFrameSecond(header);
  num=getVDIFFrameNumber(header);
  inval=getVDIFFrameInvalid(header);
 
  if(f != framebytes || inval || mjd < 50000 || mjd > 60000 || sec < 0 || num < 0 || num > 125000)
    {
      if (ifverbose)
	fprintf(stderr,"Invalid frame: bytes %i mjd %i sec %i num %i inval %i\n",f,mjd,sec,num,inval);
      return 1;
    }
  else
    return 0;
}
