#include "cvrt2to8.c"
#include "ran.c"
#include "vdif2psrfits.h"
#include <malloc.h>
#include <complex.h>
#include <fftw3.h>

void getDetection(float p0r, float p0i, float p1r, float p1i, float *det, char dstat)
{
  if(dstat=='C')
	{
	  det[0]=p0r*p0r+p0i*p0i;
	  det[1]=p1r*p1r+p1i*p1i;
	  det[2]=p0r*p1r+p0i*p1i;
	  det[3]=p0r*p1i-p0i*p1r;
	}
  else if(dstat=='I')
	det[0]=p0r*p0r+p0i*p0i+p1r*p1r+p1i*p1i;
  else if(dstat=='X')
	det[0]=p0r*p0r+p0i*p0i;
  else if(dstat=='Y')
	det[0]=p1r*p1r+p1i*p1i;
}

// Get coherence detection from real 2-bit, 32-channel frame of two pols
void getVDIFFrameDetection_32chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, float det[][4],char dstat)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts,Nchan,npol;
  unsigned char *buff8[2];

  // Decode dstat to get npol
  if(dstat=='C')
	npol=4;
  else
	npol=1;

  Nchan=32;
  for(i=0;i<2;i++)
	buff8[i]=malloc(sizeof(unsigned char)*fbytes*4);

  //Extend from 2-bit to 8-bit
  convert2to8(buff8[0], src_p0, fbytes);
  convert2to8(buff8[1], src_p1, fbytes);

  //Number of time samples
  Nts=fbytes*4/Nchan;
  
  //Memo for FFT for two pols
  in = (float **)malloc(sizeof(float *)*2);
  for(i=0;i<2;i++)
	in[i]= (float *) fftwf_malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  
  //Initialize
  for(k=0;k<Nchan;k++)
	for(j=0;j<4;j++)
	  det[k][j]=0.0;
  for(j=0;j<4;j++)
	dets[j]=0.0;
  
  //Detect each channel
  for(k=0;k<Nchan;k++)
	{
	  //For each pol
	  for(j=0;j<2;j++)
		{
		  //Read time series
		  for(i=0;i<Nts;i++)
			{
			  if(k<Nchan/2)
				in[j][i]=(float)((int)buff8[j][i*Nchan+15-k])-mean_2to8;
			  else
				in[j][i]=(float)((int)buff8[j][i*Nchan+15+Nchan-k])-mean_2to8;
			}
		}

	  //Perform FFT
	  if(k==0)
		{
		  pl0 = fftwf_plan_dft_r2c_1d(Nts, in[0], out_p0, FFTW_ESTIMATE);
		  pl1 = fftwf_plan_dft_r2c_1d(Nts, in[1], out_p1, FFTW_ESTIMATE);
		}
	  fftwf_execute(pl0);
	  fftwf_execute(pl1);
	  for(i=0;i<Nts/2+1;i++)
		{
		  getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);
		  for(j=0;j<npol;j++)
			det[k][j]+=dets[j];
		}
	}

  //Free up memo
  for(i=0;i<2;i++)
	{
	  free(buff8[i]);
	  free(in[i]);
	}
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);
}

//Generate fake detection with given mean and rms for a frame
void getVDIFFrameFakeDetection(double *mean, double *rms, int Nchan, float det[][4], long int seed, int fbytes, char dstat)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts,npol;

  // Decode dstat to get npol
  if(dstat=='C')
	npol=4;
  else
	npol=1;
  
  //Number of time samples
  Nts=fbytes*4/Nchan;
  
  //Memo for FFT for two pols
  in = (float **)malloc(sizeof(float *)*2);
  for(i=0;i<2;i++)
	in[i]= (float *)malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));

  //Initialize
  for(k=0;k<Nchan;k++)
	for(j=0;j<4;j++)
	  det[k][j]=0.0;

  for(j=0;j<4;j++)
	dets[j]=0.0;

  //Detect each channel
  for(k=0;k<Nchan;k++)
	{
	  //Operate each pol
	  for(j=0;j<2;j++)
		{
		  //Generate time series with given mean & rms
		  for(i=0;i<Nts;i++)
			  in[j][i]=mean[j]+gasdev(&seed)*rms[j];
		}
	  
	  //Perform FFT
	  if(k==0)
		{
		  pl0 = fftwf_plan_dft_r2c_1d(Nts, in[0], out_p0, FFTW_ESTIMATE);
		  pl1 = fftwf_plan_dft_r2c_1d(Nts, in[1], out_p1, FFTW_ESTIMATE);
		}
	  fftwf_execute(pl0);
	  fftwf_execute(pl1);

	  //Make detection
	  for(i=1;i<Nts/2;i++)
		{
		  getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);
		  for(j=0;j<npol;j++) 
			det[k][j]+=dets[j];
		}
	}

  //Free up memo
  for(i=0;i<2;i++)
	free(in[i]);
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);
}
