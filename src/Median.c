#include <stdio.h>
#include <math.h>

void   DSort(int , double *);

/* Simple median calculator */
double Median(double *InArray, int NPts)
{

  int i;
  double Med;
 
  DSort(NPts, InArray);
  for(i=0; i<NPts; i++)
  if(NPts%2==0){ /* even number of array elements */
    Med =  (InArray[NPts/2] + InArray[(NPts/2)-1]) / 2.;
  }
  else{  /* Odd number of array elements */
    Med =  InArray[NPts/2]; /* Will round off to the correct element */
  }

  return Med;
}

/* Runs median fileter over an array to output a smoothed array */
int MedianFilter(double *InArray, double *OutArray, int NPts, int NumNearest)
{

  int    i, j;
  double Neighbour[NPts];


  for (i=0;i<NPts;i++){
    for (j=-NumNearest/2;j<=NumNearest/2;j++){
      //printf("%d -- %d --> ",i,j);
      if(i+j < 0) { // for profile bins < half # nearest neighbours
	Neighbour[j+(NumNearest/2)] = InArray[NPts-1+j]; 
      }
      else if(i+j > NPts-1) {// for bins > (NPts - half # NN's)
	Neighbour[j+(NumNearest/2)] = InArray[j-1];
      }
      else{
	Neighbour[j+(NumNearest/2)] = InArray[i+j];
      }
      fflush(stdout);
    }
    DSort(NumNearest+1,Neighbour);  
    OutArray[i] = Neighbour[(NumNearest/2) + 1];  // choose median value
  }
  
  return 0;

}
