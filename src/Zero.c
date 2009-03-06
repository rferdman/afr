#include <math.h>
#include <stdio.h>


void IZero(int *array, int npts)
{

  int i;

  for(i=0;i<npts;i++)
    array[i] = 0;

}

void FZero(float *array, int npts)
{

  int i;

  for(i=0;i<npts;i++)
    array[i] = 0.0;

}

void DZero(double *array, int npts)
{

  int i;

  for(i=0;i<npts;i++)
    array[i] = 0.;

}

void IFill(int *array, int npts, int value)
{

  int i;

  for(i=0;i<npts;i++)
    array[i] = value;

}

void FFill(float *array, int npts, float value)
{

  int i;
  
  for(i=0;i<npts;i++)
    array[i] = value;

}
