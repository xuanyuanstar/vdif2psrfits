#include "cvrt2to8.c"
#include "ran.c"
#include <malloc.h>
#include <complex.h>
#include <fftw3.h>

void getDetection_coherence(float p0r, float p0i, float p1r, float p1i, float *det)
{
  det[0]=p0r*p0r+p0i*p0i;
  det[1]=p1r*p1r+p1i*p1i;
  det[2]=p0r*p1r+p0i*p1i;
  det[3]=p0r*p1i-p0i*p1r;
}

// Get coherence detection from real 2-bit, 32-channel frame of two pols
void getVDIFFrameDetection_coherence_32chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, float det[][4])
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts,Nchan;
  unsigned char *buff8[2];
  
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
				in[j][i]=(float)((int)buff8[j][i*Nchan+15-k]);
			  else
				in[j][i]=(float)((int)buff8[j][i*Nchan+15+Nchan-k]);
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
	  
	  //Make detection for each FFT channel and sum up
	  for(i=1;i<Nts/2+1;i++)
		{
		  //getDetection_coherence(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),-cimag(out_p1[i]),dets);
		  getDetection_coherence(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets);
		  for(j=0;j<4;j++)
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
void getFakeDetection(double *mean, double *rms, int Nchan, float det[][4], long int seed, int fbytes)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts;

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
	  for(i=1;i<Nts/2+1;i++)
		{
		  getDetection_coherence(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets);
		  for(j=0;j<4;j++) 
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

void getVDIFFrameDetection_coherence_1chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, float *det)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float **in,dets[4];
  int i,j,k,Nts;
  unsigned char *buff8[2];
  
  for(i=0;i<2;i++)
	buff8[i]=malloc(sizeof(unsigned char)*fbytes*4);

  //Extend from 2-bit to 8-bit
  convert2to8(buff8[0], src_p0, fbytes);
  convert2to8(buff8[1], src_p1, fbytes);

  //Number of time samples
  Nts=fbytes*4;

  //Memo for FFT for two pols
  in = (float **)malloc(sizeof(float *)*2);
  for(i=0;i<2;i++)
	in[i]= (float *) fftwf_malloc(sizeof(float)*Nts);
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*(Nts/2+1));
  
  //Initialize
  for(j=0;j<4;j++)
	{
	  det[j]=0.0;
	  dets[j]=0.0;
	}

  //For each pol
  for(j=0;j<2;j++)
	{
	  //Read time series
	  for(i=0;i<Nts;i++)
		in[j][i]=(float)((int)buff8[j][i]);

	  //Perform FFT
	  if(k==0)
		{
		  pl0 = fftwf_plan_dft_r2c_1d(Nts, in[0], out_p0, FFTW_ESTIMATE);
		  pl1 = fftwf_plan_dft_r2c_1d(Nts, in[1], out_p1, FFTW_ESTIMATE);
		}
	  fftwf_execute(pl0);
	  fftwf_execute(pl1);

	  //Make detection for each FFT channel and sum up
	  for(i=1;i<Nts/2+1;i++)
		{
		  //getDetection_coherence(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),-cimag(out_p1[i]),dets);
		  getDetection_coherence(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets);
		  for(j=0;j<4;j++)
			det[j]+=dets[j];
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