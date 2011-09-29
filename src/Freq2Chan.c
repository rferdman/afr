#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPDefs.h""

/* little function to map a frequency to a channel array number */
//#define TOL 0.0001


int Freq2Chan(double Freq, double *FreqArr, int n_chan) {
  
  int i_chan, Chan;
  
  for ( i_chan=0; i_chan<n_chan; i_chan++ ){
    //if(Freq == FreqArr[i_chan]) {

    /*  if(FreqArr[i_chan] < 1123. && FreqArr[i_chan] > 1121. 
       && Freq < 1123. && Freq > 1121.) 
      printf("FREQ = %lf, FREQ FREQ = %lf\nfabs = %lf\n\n", 
      Freq, FreqArr[i_chan], fabs(Freq - FreqArr[i_chan]));fflush(stdout);*/

    if(fabs(Freq - FreqArr[i_chan]) <= DBLEPS) {
      Chan = i_chan;
      return Chan;
    }
  }

  /* if we got this far, it means we didn't find the channel, so it's wrong in the file... */
  printf("Channel not found at %lf MHz.\n",Freq);
  fflush(stdout);
  return -1;

}



