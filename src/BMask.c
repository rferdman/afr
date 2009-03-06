/* Creates an array mask which denotes which bins occur 
   during the pulse (0=yes, 1=no) */ 

#include <stdio.h>
#include <math.h>
#include "ASPCommon.h"

void BMask(float *x,int *nbins,double *duty,double *y) 
{
  int    i,j,k,i0,jz,ntest,nbig;
  double ave,rms,sq,peak;

  for(i=0;i<*nbins;i++) y[i] = x[i];

  DSort(*nbins,(&y[0]-1));

  /*  for(i=0;i<(*nbins-1);i++) 
    diff[i]  =diffsave[i] = y[i+1] - y[i];
  diff[*nbins-1] = 1.e25;
  sort(*nbins,(&diff[0]-1));
  for(i=0;i<*nbins-1;i++) 
    if(diffsave[i] == diff[0]) {
      i0 = i;
      ave = y[i];
      printf("Ave, i0: %f, %d\n",ave,i0);
      break;
    } */
  i0 = (0.5*(*duty)*(*nbins)+0.5);
  ave= y[i0-1];  
  peak= y[*nbins-1];

  sq = 0.0;
  for(i=0;i<i0;i++) sq += (y[i] - ave) * (y[i] - ave);
  rms = sqrt(sq/i0);

  if((*nbins/128) > 2) jz = *nbins/128;
  else jz = 2;
  ntest = jz/2 - 1;

  for(i=0;i<*nbins;i++) {
    y[i] = 1.;
    nbig = 0;
    for(j=-jz;j<jz+1;j++) {
      k = i + j;
      if(k<0) k += *nbins;
      if(k>*nbins) k -= *nbins;
      if((x[k] - ave) > 1.8*rms) nbig++;
    }
    if(nbig > ntest) y[i] = 0.;
  }
}


/* Calculates the rms signal for the baseline */
void Baseline(float *x,double *y, int *nbins, double *ave, double *rms)
{

  double s,sq,sw;
  int i;

  s = sw = sq = 0.;

  for(i=0;i<*nbins;i++) {
    s += x[i]*y[i];
    sw += y[i];
  }
  *ave = s/sw;
  /*  printf("sw = %4.0f\n",sw);  */

  for(i=0;i<*nbins;i++)
    sq += y[i] * (x[i] - *ave) * (x[i] - *ave);
  *rms = sqrt(sw/sq);   /* Want 1/rms  */

}


