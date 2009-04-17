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

  for (i_arr=0; i_arr<ArrLen; i_arr++)
    /* As soon as we find one bin not zeroed, return with ArrayZero=0 */
    if (Array[i_arr]!=0.0) return 0;

  /* If we made it this far, then all bins have value of 0.0, so ArrayZero=1 */
  return 1 ;
  
  

}


