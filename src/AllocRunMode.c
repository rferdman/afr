/* Will dnamically allocate all non-char array members of the RunMode 
   structure, so as to avoid memory leaks, maximum array sizes, etc.

   Robert Ferdman, 20 April 2010 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPDefs.h"
#include "ASPCommon.h"

int AllocRunMode (struct RunVars *RunMode) {

  RunMode->FirstChanAdd = (int *) malloc (NCHMAX*sizeof(int));
  RunMode->LastChanAdd =  (int *) malloc (NCHMAX*sizeof(int));
  RunMode->MinChans2Add = (int *) malloc (NCHMAX*sizeof(int));
  RunMode->MaxChans2Add = (int *) malloc (NCHMAX*sizeof(int));
  RunMode->NumChans2Add = (int *) malloc (NCHMAX*sizeof(int));
  RunMode->CurZapChan =   (int *) malloc (NCHMAX*sizeof(int));
  RunMode->ZapChan =      (int *) malloc (NCHMAX*sizeof(int));
  RunMode->DumpOmit =     (int *) malloc (MAXOMIT*sizeof(int));
  RunMode->ChanOmit =     (int *) malloc (MAXOMIT*sizeof(int));
  RunMode->TotScans =     (int *) malloc (MAXDUMPS*NCHMAX*sizeof(int));
  RunMode->TotOmit =      (int *) malloc (MAXDUMPS*NCHMAX*sizeof(int));
  RunMode->AllDumpOmit =  (int *) malloc (MAXOMIT*sizeof(int));
  RunMode->OmitFlag =     (int *) malloc (MAXDUMPS*NCHMAX*sizeof(int));
  RunMode->MM =         (float *) malloc (16*NCHMAX*sizeof(float));
  RunMode->ThetaBB =   (double *) malloc (NCHMAX*sizeof(double));


  return 0;

}
