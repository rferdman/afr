#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* little function to tell whether all elements of an array have a value
   of zero.  Returns 0 if even one element is non-zero, and 1 if all have 
   been zeroed out */


int ArrayZero(float *Array, int ArrLen)
{
  int i_arr;
  int not_all_zeroes=0;

  for (i_arr=0; i_arr<ArrLen; i_arr++){
    
    /* First check if the array element is a fintite, or is NaN.  
       If so, just return a 2 saying that this array is no good */
    if (isnan(Array[i_arr])) return 2;
    if (!isfinite(Array[i_arr])) return 3;
    /* As soon as we find one bin not zeroed, Flick the not_all_zeroes switch */
    if (Array[i_arr]!=0.0) not_all_zeroes=1;
  }

  /* If we made it this far, there are no nans or infiinite numbers */
  /* Retuen 0 if the array is not all zeroes, and 1 if it is all zeroes */
  if(not_all_zeroes) 
    return 0;
  else 
    return 1 ;
  
  

}


