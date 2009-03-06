/* Function that reads in a set of Stokes profile and adds user-defined dumps together
   cumulatively.  Since we do not need to rotate the profiles, it is a simple case of 
   adding arrays together.*/

/* Output profiles are the added profiles, as well as the non-added profiles */

/* --RDF 21/04/04 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


void AddScans(struct RunVars *RunMode, struct ASPHdr *hdr, 
	      struct SubHdr *SubInHdr, struct SubHdr *SubDumpHdr,
	      struct StdProfs *InputProfs, struct StdProfs *OutputProfs)
{

  
  int             i,j,k, nomit, totalzapchan;
  int             SubHdrIndex, curomit[NCHMAX];
  static int      NumDumps=0, omitcount[NCHMAX];
  struct StdProfs TempProfs, WriteProfs;
  double          denom;
  char            Outputfile[256];
  
  








  for (i_dump=0; i_dump<RunMode->NOutDumps; i_dump++){
    /* If we are in the last Output dump, max dump to add is actual last 
       dump so we don't go over and add things that don't exist */
    MinAddDump = i_dump*RunMode->AddDumps; 
    if(i_dump == RunMode->NOutDumps-1){
      MaxAddDump = RunMode->NDumps; // not NDumps-1 to include last dump
    }
    else{
      //	MaxAddDump = i_dump*RunMode->AddDumps + RunMode->AddDumps;
      MaxAddDump = MinAddDump + RunMode->AddDumps;
    }





    
    
    Zero(ThisDumpProf);
    
    /* Do one (output) channel at a time */
    for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
           







      Zero(NewProf);

      /* Now loop over beginning and end dumps of each output dump */
      for (i_dump_add=MinAddDump; i_dump_add<MaxAddDump; i_dump_add++){
       	
	








	/* Now loop over beginning and end channels of each output channel */
	/* Note that we put "<=" here since we want the range to be 
	   inclusive */
	for (i_chan_add=RunMode->FirstChanAdd[i_chan]; 
	     i_chan_add<=RunMode->LastChanAdd[i_chan]; i_chan_add++){ 
	  
	  NewProf += OldProf[i_dump_add][i_chan_add];








	}










      }
      








      ThisDumpProf[i_chan] += NewProf;
    }












    WrtASPStokes(ThisDumpProf);


    free(NewProf);
    free(ThisDumpProf);

  }



	  RunMode->TotScans[i_dump*RunMode->NOutChans + i_chan]++;
	  RunMode->TotOmit[i_dump*RunMode->NOutChans + i_chan] += 
	    RunMode->OmitFlag[i_dump_add*hdr->obs.NChan + i_chan_add];









  /*************************** OLD STUFF BEGINS ****************************/

  /* If not adding dumps, then just straight copy */
  if(RunMode->AddDumps == 1){

    /* Check to see if nscan is equal to a dump we want to omit entirely */
    for (k=0;k<RunMode->NumAllDumpOmit;k++){
      if (nscan == RunMode->AllDumpOmit[k]){
	*NewDump = 0;
	printf("Dump %d is omitted entirely...\n",nscan);fflush(stdout);
	return;
      }
    }

    SubDumpHdr->DumpMiddleSecs = SubInHdr[nscan].DumpMiddleSecs;
    for(i=0;i<hdr->obs.NChan;i++){
      SubDumpHdr->DumpRefPhase[i] = SubInHdr[nscan].DumpRefPhase[i];
      SubDumpHdr->DumpRefPeriod[i] = SubInHdr[nscan].DumpRefPeriod[i];
      memcpy(&OutputProfs[i], &InputProfs[i], sizeof(struct StdProfs));
    }
    *NewDump = 1; // make zero if you don't want to include entire dump in file?
  }
  else{  
    
    if((nscan == 0) || (nscan%RunMode->AddDumps == 0)){
      /* First time through or start of new dump, 
       * just zero out profiles to start */
      for(i=0;i<hdr->obs.NChan;i++){
	FZero(OutputProfs[i].rstds,NBINMAX); 
	FZero(OutputProfs[i].rstdq,NBINMAX);
	FZero(OutputProfs[i].rstdu,NBINMAX); 
	FZero(OutputProfs[i].rstdv,NBINMAX); 
	FZero(OutputProfs[i].stdlin,NBINMAX); 
	FZero(OutputProfs[i].stdphi,NBINMAX); 
	FZero(OutputProfs[i].stdphierr,NBINMAX);
	omitcount[i]=0;
      }
    }

    /* accumulate profile */
    for(i=0;i<hdr->obs.NChan;i++){
	
      /* See if any omissions to take into account */
      curomit[i]=0;
      for(nomit=0;nomit<RunMode->NScanOmit;nomit++){
	if((nscan==RunMode->DumpOmit[nomit]) && (i==RunMode->ChanOmit[nomit])){
	  curomit[i]=1;
	  omitcount[i]++;
	}
      }
      for (k=0;k<RunMode->NumAllDumpOmit;k++){
	if (nscan == RunMode->AllDumpOmit[k]){     
	  curomit[i]=1;
	  omitcount[i]++;
	  if (i==0)
	    printf("Dump %d is omitted entirely...\n",nscan);fflush(stdout);
	  break;
	}
      }
      
      if(!curomit[i]){
	for(j=0;j<RunMode->NBins;j++) {
	  OutputProfs[i].rstds[j] += InputProfs[i].rstds[j];
	  OutputProfs[i].rstdq[j] += InputProfs[i].rstdq[j];
	  OutputProfs[i].rstdu[j] += InputProfs[i].rstdu[j];
	  OutputProfs[i].rstdv[j] += InputProfs[i].rstdv[j];
	}
      }

    }

    /* If it is dump time... */
    if((nscan>0) 
       && (nscan != RunMode->NDumps - 1) 
       && ((nscan+1)%RunMode->AddDumps == 0)) {

      if(RunMode->Verbose)
	printf("putting together dump number %d\n",NumDumps);fflush(stdout);

    /*  if (nscan == RunMode->AllDumpOmit[k]){
        *NewDump = 0;
        printf("This dump (would be %d) is omitted entirely due to\n",NumDumps);
        printf("lack of scans going into it...\n");fflush(stdout);
        NumDumps++;
        return;
      }  */

      SubHdrIndex = (int)(floor((double)((NumDumps*RunMode->AddDumps + 
			       ((NumDumps+1)*(RunMode->AddDumps) - 1))/2.)));
      SubDumpHdr->DumpMiddleSecs = SubInHdr[SubHdrIndex].DumpMiddleSecs;

      totalzapchan = 0;
      for(i=0;i<hdr->obs.NChan;i++){

        RunMode->CurZapChan[i] = 0;
        if((RunMode->AddDumps - omitcount[i]) <= 0) {
          RunMode->CurZapChan[i] = 1;
          totalzapchan++; 
          printf("Bad denominator excluded from dump %d, channel %d\n",
                  NumDumps,i);fflush(stdout);
        }

	sprintf(Outputfile,"%s.%4.4d.%4.4d.out.dump",
		RunMode->OutfileRoot,(int)(hdr->obs.ChanFreq[i]),NumDumps);
	if (RunMode->Verbose)
	  printf("Outputfile = %s\n",Outputfile);fflush(stdout);

	SubDumpHdr->DumpRefPhase[i] = SubInHdr[SubHdrIndex].DumpRefPhase[i];
	SubDumpHdr->DumpRefPeriod[i] = SubInHdr[SubHdrIndex].DumpRefPeriod[i];

	/* Average over number of dumps added together, taking into account
	 * whether any dump/freq combos were omitted*/
	for(k=0;k<RunMode->NBins;k++){
	  TempProfs.rstds[k]=OutputProfs[i].rstds[k]/
	    ((double)(RunMode->AddDumps - omitcount[i]));
	  TempProfs.rstdq[k]=OutputProfs[i].rstdq[k]/
	    ((double)(RunMode->AddDumps - omitcount[i]));
	  TempProfs.rstdu[k]=OutputProfs[i].rstdu[k]/
	    ((double)(RunMode->AddDumps - omitcount[i]));
	  TempProfs.rstdv[k]=OutputProfs[i].rstdv[k]/
	    ((double)(RunMode->AddDumps - omitcount[i]));
	}
	if(RunMode->BinDown) {
	  BinDown(RunMode, &TempProfs, &WriteProfs); 
	}
	else {
	  memcpy(&WriteProfs, &TempProfs, sizeof(struct StdProfs));
	}
	MakePol(RunMode, RunMode->NBinsOut, &WriteProfs);

	memcpy(&OutputProfs[i],&TempProfs, sizeof(struct StdProfs));
	OutputProfs[i].SNR = WriteProfs.SNR;

	/* WriteStokes(RunMode, &WriteProfs, OutputHead[i], Outputfile); */
      }
      if(totalzapchan == hdr->obs.NChan){ 
        *NewDump = 0;
        printf("All channels in dump %d have bad denom\n",
               NumDumps);fflush(stdout);
      }
      else
        *NewDump = 1;
        
      NumDumps++;


    }
    
    
    
  
    //  } //put this under finaldump???

  /* Put in last dump here */
  if(nscan == RunMode->NDumps - 1){
    if (RunMode->NDumps%RunMode->AddDumps == 0)
      denom = RunMode->AddDumps;
    else
      denom = RunMode->NDumps%RunMode->AddDumps;

    SubHdrIndex = (int)(floor((double)((RunMode->NumEffDumps-1)*RunMode->AddDumps
				       + (RunMode->NDumps-1))/2.));
    SubDumpHdr->DumpMiddleSecs = SubInHdr[SubHdrIndex].DumpMiddleSecs;

    totalzapchan = 0;
    for(i=0;i<hdr->obs.NChan;i++){

      RunMode->CurZapChan[i] = 0;
      if((denom - omitcount[i]) <= 0) {
        RunMode->CurZapChan[i] = 1;
        totalzapchan++; 
      }

      sprintf(Outputfile,"%s.%4.4d.%4.4d.out.dump",
	      RunMode->OutfileRoot,(int)(hdr->obs.ChanFreq[i]),
	      RunMode->NumEffDumps-1);
      if(RunMode->Verbose)
	printf("Outputfile = %s\n",Outputfile);fflush(stdout);

      SubDumpHdr->DumpRefPhase[i] = SubInHdr[SubHdrIndex].DumpRefPhase[i];
      SubDumpHdr->DumpRefPeriod[i] = SubInHdr[SubHdrIndex].DumpRefPeriod[i];
      
      /* Average over number of dumps added together, taking into account
       * whether any dump/freq combos were omitted */
      for(k=0;k<RunMode->NBins;k++){
	TempProfs.rstds[k] = OutputProfs[i].rstds[k] /
	  (double)(denom - omitcount[i]);
	TempProfs.rstdq[k] = OutputProfs[i].rstdq[k] /
	  (double)(denom - omitcount[i]);
	TempProfs.rstdu[k] = OutputProfs[i].rstdu[k] /
	  (double)(denom - omitcount[i]);
	TempProfs.rstdv[k] = OutputProfs[i].rstdv[k] /
	  (double)(denom - omitcount[i]);
      }
      memcpy(&OutputProfs[i], &TempProfs,sizeof(struct StdProfs));

      if(RunMode->BinDown){
	BinDown(RunMode, &TempProfs, &WriteProfs);
      }
      else {
	memcpy(&WriteProfs, &TempProfs, sizeof(struct StdProfs));
      }
      MakePol(RunMode, RunMode->NBinsOut, &WriteProfs);

      OutputProfs[i].SNR = WriteProfs.SNR;
	  
      /*  WriteStokes(RunMode, &WriteProfs, OutputHead[i], Outputfile); */
    }
    if(totalzapchan == hdr->obs.NChan) 
      *NewDump = 0;
    else
      *NewDump = 1;

  }
  
  } //put this under finaldump???


}
