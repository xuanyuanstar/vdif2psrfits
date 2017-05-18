//Calculate totals of odd and even samples in a 2-bit VDIF frame, for odd and even samples
#include "cvrt2to8.c"
#include <malloc.h>
#include <stdio.h>

// Single channel
void getVDIFFrameTotal(const unsigned char *src, int fbytes, int *tot)
{
  int i;
  unsigned char *buff8;

  buff8=malloc(sizeof(unsigned char)*fbytes*4);
  tot[0]=0;
  tot[1]=0;

  //Extend from 2-bit to 8-bit
  convert2to8(buff8, src, fbytes);

  for(i=0;i<fbytes*2;i++)
	{
	  tot[0]+=(int)buff8[2*i];
	  tot[1]+=(int)buff8[2*i+1];
	}
  free(buff8);
}

// 32 channels, 2-bit real sampled
void getVDIFFrameTotal_32chan(const unsigned char *src, int fbytes, int tot[][32])
{
  int i,j,k,Ncyc,Nchan,**dat;
  unsigned char *buff8;

  Nchan=32;
  
  dat=(int **)malloc(sizeof(int *)*2);
  dat[0]=(int *)malloc(sizeof(int)*Nchan);
  dat[1]=(int *)malloc(sizeof(int)*Nchan);
  
  buff8=malloc(sizeof(unsigned char)*fbytes*4);

  // Initialize
  for(i=0;i<2;i++)
	for(j=0;j<Nchan;j++)
		tot[i][j]=0;

  //Extend from 2-bit to 8-bit
  convert2to8(buff8, src, fbytes);

  //Number of cycles of two time samples
  Ncyc=fbytes*4/Nchan/2;

  for(i=0;i<Ncyc;i++)
	{
	  for(j=0;j<2;j++)
		{
		  for(k=0;k<Nchan;k++)
			{
			  if(k<16)
				{
				  tot[j][k]+=(int)buff8[i*Nchan*2+j*Nchan+15-k];
				}
			  else				
				{
				  tot[j][k]+=(int)buff8[i*Nchan*2+j*Nchan+15-k+Nchan];
				}
			}
		}
	}
  free(buff8);
}

void getDetection_coherence(int p0r, int p0i, int p1r, int p1i, int *det)
{
  det[0]=p0r*p0r+p0i*p0i;
  det[1]=p1r*p1r+p1i*p1i;
  det[2]=p0r*p1r+p0i*p1i;
  det[3]=p0r*p1i-p0i*p1r;
}
  

void getVDIFFrameDetection_coherence_32chan(const unsigned char *src_p0, const unsigned char *src_p1, int fbytes, int det[][4])
{
  int i,j,k,Ncyc,Nchan,det_s[4];
  unsigned char *buff8[2];

  Nchan=32;
  
  for(i=0;i<2;i++)
	buff8[i]=malloc(sizeof(unsigned char)*fbytes*4);
  
  //Extend from 2-bit to 8-bit
  convert2to8(buff8[0], src_p0, fbytes);
  convert2to8(buff8[1], src_p1, fbytes);

  //Number of cycles of two time samples   
  Ncyc=fbytes*4/Nchan/2;

  for(k=0;k<Nchan;k++)
	for(j=0;j<4;j++)
	  det[k][j]=0;
  
  for(i=0;i<Ncyc;i++)
	{
	  for(k=0;k<Nchan;k++)
		{
		  if(k<16)
			{
			  getDetection_coherence((int)buff8[0][i*Nchan*2+15-k],(int)buff8[0][i*Nchan*2+Nchan+15-k],(int)buff8[1][i*Nchan*2+15-k],(int)buff8[1][i*Nchan*2+Nchan+15-k],det_s);
			  for(j=0;j<4;j++)
				det[k][j]+=det_s[j];
				  //det[k][0]=(int)buff8[0][i*Nchan*2+15-k]*(int)buff8[0][i*Nchan*2+15-k]+(int)buff8[0][i*Nchan*2+Nchan+15-k]*(int)buff8[0][i*Nchan*2+Nchan+15-k];
				  //det[k][1]=(int)buff8[1][i*Nchan*2+15-k]*(int)buff8[1][i*Nchan*2+15-k]+(int)buff8[1][i*Nchan*2+Nchan+15-k]*(int)buff8[1][i*Nchan*2+Nchan+15-k];
				  //det[k][2]=(int)buff8[0][i*Nchan*2+15-k]*(int)buff8[1][i*Nchan*2+15-k]+(int)buff8[0][i*Nchan*2+Nchan+15-k]*(int)buff8[0][i*Nchan*2+Nchan+15-k];
				  //det[k][3]=(int)buff8[0][i*Nchan*2+15-k]*(int)buff8[1][i*Nchan*2+Nchan+15-k]-(int)buff8[0][i*Nchan*2+Nchan+15-k]*(int)buff8[1][i*Nchan*2+15-k];
				  //(int)buff8[0][i*Nchan*2+15-k];
				  //(int)buff8[0][i*Nchan*2+Nchan+15-k];
				  //(int)buff8[1][i*Nchan*2+15-k];
				  //(int)buff8[1][i*Nchan*2+Nchan+15-k];				  
			}
		  else
			{
			  getDetection_coherence((int)buff8[0][i*Nchan*2+15-k+Nchan],(int)buff8[0][i*Nchan*2+Nchan+15-k+Nchan],(int)buff8[1][i*Nchan*2+15-k+Nchan],(int)buff8[1][i*Nchan*2+Nchan+15-k+Nchan],det_s);
			  for(j=0;j<4;j++)
				det[k][j]+=det_s[j];
				  //(int)buff8[0][i*Nchan*2+15-k+Nchan];
				  //(int)buff8[0][i*Nchan*2+Nchan+15-k+Nchan];
				  //(int)buff8[1][i*Nchan*2+15-k+Nchan];
				  //(int)buff8[1][i*Nchan*2+Nchan+15-k+Nchan];				  
			}
		}
	}
  //for(k=0;k<Nchan;k++)
  //for(j=0;j<4;j++)
  //printf("%i\n",det[k][j]);
  
  for(i=0;i<2;i++)
	free(buff8[i]);
}
