#include <stdio.h>
#include <math.h>
                                                                                
                                                                                
double Max(double *array,int nbins, int *imax)
{
                                                                                
  int    i;
  double max;
                                                                                
  max = -999999.;
  *imax = 0;
                                                                                
  for(i=0;i<nbins;i++)
    if(array[i] > max) {
      max = array[i];
      *imax = i;
    }
                                                                                
  return(max);
                                                                                
}

double Min(double *array,int nbins, int *imin)
{
                                                                                
  int    i;
  double min;
                                                                                
  min = array[0];
  *imin = 0;
                                                                                
  for(i=1;i<nbins;i++)
    if(array[i] < min) {
      min = array[i];
      *imin = i;
    }
                                                                                
  return(min);
                                                                                
}

                                                                                
float FMax(float *array,int nbins, int *imax)
{
                                                                                
  int    i;
  float max;
                                                                                
  max = -999999.;
  *imax = 0;
                                                                                
  for(i=0;i<nbins;i++)
    if(array[i] > max) {
      max = array[i];
      *imax = i;
    }
                                                                                
  return(max);
                                                                                
}

float FMin(float *array,int nbins, int *imin)
{
                                                                                
  int    i;
  float min;
                                                                                
  min = array[0];
  *imin = 0;
                                                                                
  for(i=1;i<nbins;i++)
    if(array[i] < min) {
      min = array[i];
      *imin = i;
    }
                                                                                
  return(min);
                                                                                
}

int IMax(int *array,int nbins, int *imax)
{
                                                                                
  int    i;
  int    max;
                                                                                
  max = -99999;
  *imax = 0;
                                                                                
  for(i=0;i<nbins;i++)
    if(array[i] > max) {
      max = array[i];
      *imax = i;
    }
                                                                                
  return(max);
                                                                                
}

int IMin(int *array,int nbins, int *imin)
{
                                                                                
  int    i;
  int    min;
                                                                                
  min = array[0];
  *imin = 0;
                                                                                
  for(i=1;i<nbins;i++)
    if(array[i] < min) {
      min = array[i];
      *imin = i;
    }
                                                                                
  return(min);
                                                                                
}

