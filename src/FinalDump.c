/* Function that does the final dump, if NDumps%(Number of dumps to add together 
   at a time !=0).  Since we do not need to rotate the profiles, it is a simple 
   case of adding arrays together.*/

/* Output profiles are the added profiles, as well as the non-added profiles */

/* --RDF 23/05/04 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


void FinalDump(struct RunVars *RunMode, struct ASPHdr *hdr, 
	       struct StdProfs *OutputProfs, char **OutputHead)
{

 
  int i,k;
  double denom;
  struct StdProfs TempProfs, WriteProfs;
  char Outputfile[256];

  /* HERE: LOOP over number of sets of additions to make */  
  //  printf("RunMode->AddDumps = %d\n",RunMode->AddDumps);fflush(stdout);

  if (RunMode->NDumps%RunMode->AddDumps == 0)
    denom = (double)(RunMode->AddDumps);
  else
    denom = (double)(RunMode->NDumps%RunMode->AddDumps);

  for(i=0;i<hdr->obs.NChan;i++){


    sprintf(Outputfile,"%s.%4.4d.%4.4d.out.dump",
	    RunMode->OutfileRoot,(int)(hdr->obs.ChanFreq[i]),
	    RunMode->NumEffDumps-1);
    if(RunMode->Verbose)
      printf("Outputfile = %s\n",Outputfile);fflush(stdout);

    for(k=0;k<RunMode->NBins;k++){
      TempProfs.rstds[k] = OutputProfs[(RunMode->NumEffDumps-1)*
				       hdr->obs.NChan + i].rstds[k]/denom;
      TempProfs.rstdq[k] = OutputProfs[(RunMode->NumEffDumps-1)*
				       hdr->obs.NChan + i].rstdq[k]/denom;
      TempProfs.rstdu[k] = OutputProfs[(RunMode->NumEffDumps-1)*
				       hdr->obs.NChan + i].rstdu[k]/denom;
      TempProfs.rstdv[k] = OutputProfs[(RunMode->NumEffDumps-1)*
				       hdr->obs.NChan + i].rstdv[k]/denom;
    }
    memcpy(&OutputProfs[(RunMode->NumEffDumps-1)*hdr->obs.NChan + i],
	   &TempProfs,sizeof(struct StdProfs));

    if(RunMode->BinDown){
      BinDown(RunMode, &TempProfs, &WriteProfs);
    }
    else {
      memcpy(&WriteProfs, &TempProfs, sizeof(struct StdProfs));
    }
    MakePol(RunMode, RunMode->NBinsOut, &WriteProfs);

    OutputProfs[(RunMode->NumEffDumps-1)*hdr->obs.NChan + i].SNR =WriteProfs.SNR;
	  
   /*  WriteStokes(RunMode, &WriteProfs, OutputHead[i], Outputfile); */
  


	  


  }
      
}
