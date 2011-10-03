#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPDefs.h""

/* little function to map a frequency to a channel array number */

int Freq2Chan(double Freq, double *FreqArr, int n_chan) {
  
  int i_chan, Chan;
  
  for ( i_chan=0; i_chan<n_chan; i_chan++ ){


    if(fabs(Freq - FreqArr[i_chan]) <= DBLEPS) {
      Chan = i_chan;
      return Chan;
    }
  }

  /* if we got this far, it means we did not find the channel, so it is wrong in the file... */
  printf("Channel not found at %lf MHz.\n",Freq);
  fflush(stdout);
  return -1;

}



