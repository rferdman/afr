#include <stdio.h>
#include <math.h>


double FindPeak(float *prof,int *nbins, int *ipk)
{

  int    i;
  double max;

  max = -999999.;
  *ipk = 0;

  for(i=0;i<*nbins;i++) 
    if(prof[i] > max) {
      max = prof[i];
      *ipk = i;
    }

  return(max);

}
