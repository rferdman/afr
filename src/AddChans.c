/* Function that reads in a set of Stokes profile and adds user-defined channels
   together.  It does this for a given set of channels by rotating every channel 
   to match the phase of the reference profile's fundamental. */

/* Output profiles are the added profiles, as well as the non-added profiles */

/* --RDF 10/03/04 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

void AddChans(struct RunVars *RunMode, struct StdProfs *InputProfs, 
	      struct StdProfs *OutputProfs, struct StdProfs *StdProfile, 
	      char **HeadLine, struct ASPHdr *hdr,
	      struct SubHdr *SubDumpHdr, struct SubHdr *SubOutHdr,
	      double *OutputFreq, char **OutputHead, int nscan)
{

  int             i,j,k, zapchans, nomit;
  int             curomit, omitcount;
  int             ChanIndex, RefIndex[NCHMAX], notpicked;
  double          ByPhase, RefFreq, RefPeriod;
  struct StdProfs *AddProfs, *TempProfs;

  int test_file=1, i_chan;
  char testfile[32];
  FILE *Ftest;


  /* Malloc'ing */
  AddProfs = (struct StdProfs *)malloc(NCHMAX*NDSMAX*sizeof(struct StdProfs));  
  TempProfs = (struct StdProfs *)malloc(hdr->obs.NChan*sizeof(struct StdProfs));  


  /* LOOP over number of sets of additions to make */  

  if(RunMode->AddChans){
    for(i=0;i<RunMode->SetsofChans2Add;i++){
      zapchans=0;

      /* Zero the added profile */
      FZero(AddProfs[i].rstds,NBINMAX);
      FZero(AddProfs[i].rstdq,NBINMAX);
      FZero(AddProfs[i].rstdu,NBINMAX);
      FZero(AddProfs[i].rstdv,NBINMAX);
      FZero(AddProfs[i].stdlin,NBINMAX);
      FZero(AddProfs[i].stdphi,NBINMAX);
      FZero(AddProfs[i].stdphierr,NBINMAX);


      /* Find the phase, from reference, TO WHICH all others must rotate */
      /* Make it the middle channel out of the ones added */
      RefIndex[i] = (int)(((double)(RunMode->MinChans2Add[i] + 
				    RunMode->MaxChans2Add[i]))/2.);
      

      cprofc(InputProfs[RefIndex[i]].rstds,RunMode->NBins,
	     InputProfs[RefIndex[i]].stdamp,InputProfs[RefIndex[i]].stdpha); 


      //      RefFreq = hdr->obs.ChanFreq[RefIndex[i]]; 
      //      RefPeriod = SubDumpHdr->DumpRefPeriod[RefIndex[i]];


      omitcount=0;
      /* Rotate each profile in set to match this reference phase */
      for(j=RunMode->MinChans2Add[i];j<RunMode->MaxChans2Add[i]+1;j++) {

	curomit=0;
	/* Omit scans only if we are NOT adding dumps */
	if (RunMode->AddDumps == 1){
	  /* Cycle through each scan/channel combination to be omitted */
	  for(nomit=0;nomit<RunMode->NScanOmit;nomit++){
	    /* if the current scan and channel combo are to be omitted, flag */
	    if((nscan==RunMode->DumpOmit[nomit]) && 
	       (j==RunMode->ChanOmit[nomit])){
	      curomit=1;
	      omitcount++;
	    }
	  }
	}
	
	if (!curomit) {	  /* if scan/chan combo is NOT to be omitted */

	  if (j != RefIndex[i]){ 

	    // can probably get rid of this empty if statement 
	    
	    /* Write out individual rotated profiles if asked for */
	  
	  }

	  if(!RunMode->ZapChan[j]){
          if(!RunMode->CurZapChan[j]){

	    /* HERE: add to profile */
	    for(k=0;k<RunMode->NBins;k++){
	      AddProfs[i].rstds[k] += InputProfs[j].rstds[k];
	      AddProfs[i].rstdq[k] += InputProfs[j].rstdq[k];
	      AddProfs[i].rstdu[k] += InputProfs[j].rstdu[k];
	      AddProfs[i].rstdv[k] += InputProfs[j].rstdv[k];
	    }
	  
	  }
          }
	  else{
	    zapchans++;
	  }

	} 
	else{
	  //omitcount++;  
	}
	
      }

      /* HERE: average over number of channels added to make profile */
      /*       *** eventually weighted average? ***                  */
      for(k=0;k<RunMode->NBins;k++){
	AddProfs[i].rstds[k] /= (RunMode->NumChans2Add[i] - zapchans - omitcount);
	AddProfs[i].rstdq[k] /= (RunMode->NumChans2Add[i] - zapchans - omitcount);
	AddProfs[i].rstdu[k] /= (RunMode->NumChans2Add[i] - zapchans - omitcount);
	AddProfs[i].rstdv[k] /= (RunMode->NumChans2Add[i] - zapchans - omitcount);
      }

      // print out denominator here...
      if (test_file){
	sprintf(testfile,"test_denom.dat");
	Ftest=fopen(testfile, "a");
	fprintf(Ftest, "ADDCHANS DENOMINATOR -- DUMP %d:\n\n",nscan);
      fprintf(Ftest,"Combined Channel %d: \n", i);
      fprintf(Ftest," %d      --->  %d\n",
	      RunMode->MinChans2Add[i],RunMode->MaxChans2Add[i]);
      fprintf(Ftest," %6.1lf  --->  %6.1lf\n", 
	     hdr->obs.ChanFreq[RunMode->MinChans2Add[i]], 
	     hdr->obs.ChanFreq[RunMode->MaxChans2Add[i]]);
	fprintf(Ftest, "NumChans2Add = %d\n",RunMode->NumChans2Add[i]);
	fprintf(Ftest, "zapchans = %d\n",zapchans);
	fprintf(Ftest, "omitcount = %d\n",omitcount);
	fprintf(Ftest, "DENOMINATOR = %d\n\n\n",RunMode->NumChans2Add[i] - zapchans - omitcount);
	fclose(Ftest);
      }


    }

    /* Now i have RunMode->SetsofChans2Add number of added-together profiles */

    /* Put into Output Profs, with correct Headline, including non-added profs */
    /* and bin down (if desired), create stdlin, stdphi, stdphierr  */
    RunMode->NumEffChans = 0;  /* Initialize Number of "effective" channels */
 
    if(RunMode->Verbose)
      printf("\nCHANNELS TO ADD:\n");fflush(stdout);

    /* Pick out channels not combined: */
    /* Use min channel index as max for initial "not picked" channel */
    notpicked = RunMode->MinChans2Add[0]; 

    for(j=0;j<notpicked;j++){
      if(RunMode->Verbose)
	printf("Channel %d  OR  %f MHz --> going solo\n",
	       j,hdr->obs.ChanFreq[j]);fflush(stdout);
      
      ChanIndex = j;  /* Determine original channel index */

      /* Assign Output Profile to each individual channel */
      for(k=0;k<RunMode->NBins;k++){
	TempProfs[RunMode->NumEffChans].rstds[k] =InputProfs[ChanIndex].rstds[k];
 	TempProfs[RunMode->NumEffChans].rstdq[k] =InputProfs[ChanIndex].rstdq[k];
	TempProfs[RunMode->NumEffChans].rstdu[k] =InputProfs[ChanIndex].rstdu[k];
	TempProfs[RunMode->NumEffChans].rstdv[k] =InputProfs[ChanIndex].rstdv[k];
      }
      /* Assign central frequencies */
      OutputFreq[RunMode->NumEffChans] = hdr->obs.ChanFreq[ChanIndex];
      OutputHead[RunMode->NumEffChans] = HeadLine[ChanIndex];
      SubOutHdr->DumpRefPhase[RunMode->NumEffChans] = 
	SubDumpHdr->DumpRefPhase[ChanIndex];
      SubOutHdr->DumpRefPeriod[RunMode->NumEffChans] =
	SubDumpHdr->DumpRefPeriod[ChanIndex];
      RunMode->NumEffChans++;  /* Increase number of effective channels */

    }

    for(i=0;i<RunMode->SetsofChans2Add;i++){
      if(RunMode->Verbose)
	printf("Channels %d to %d  OR  %f to %f MHz\n",
	       RunMode->MinChans2Add[i],RunMode->MaxChans2Add[i],
	       hdr->obs.ChanFreq[RunMode->MinChans2Add[i]],
	       hdr->obs.ChanFreq[RunMode->MaxChans2Add[i]]);fflush(stdout);


      /* Assign Output Profile to added-together channel */
      for(k=0;k<RunMode->NBins;k++){
	TempProfs[RunMode->NumEffChans].rstds[k] = AddProfs[i].rstds[k];
 	TempProfs[RunMode->NumEffChans].rstdq[k] = AddProfs[i].rstdq[k];
	TempProfs[RunMode->NumEffChans].rstdu[k] = AddProfs[i].rstdu[k];
	TempProfs[RunMode->NumEffChans].rstdv[k] = AddProfs[i].rstdv[k];
      }
      /* Assign central frequency */
      OutputFreq[RunMode->NumEffChans] = hdr->obs.ChanFreq[RefIndex[i]];
      OutputHead[RunMode->NumEffChans] = HeadLine[RefIndex[i]];
      SubOutHdr->DumpRefPhase[RunMode->NumEffChans] = 
	SubDumpHdr->DumpRefPhase[RefIndex[i]];
      SubOutHdr->DumpRefPeriod[RunMode->NumEffChans] =
	SubDumpHdr->DumpRefPeriod[RefIndex[i]];
      if(RunMode->Verbose)
	printf("   -- Effective channel is %d  OR  %f MHz\n",RunMode->NumEffChans,
	       OutputFreq[RunMode->NumEffChans]);fflush(stdout);
      RunMode->NumEffChans++;  /* Increase number of effective channels */
	
      /* Pick out channels not combined: */
      if(i<RunMode->SetsofChans2Add-1){
	notpicked = RunMode->MinChans2Add[i+1] - RunMode->MaxChans2Add[i] - 1;

  
	for(j=0;j<notpicked;j++){
	  if(RunMode->Verbose){
	    printf("Channel %d  OR  %f MHz --> going solo\n",
		   j+RunMode->MaxChans2Add[i]+1, 
		   hdr->obs.ChanFreq[j+RunMode->MaxChans2Add[i]+1]);
	    fflush(stdout);
	  }

	  /* Determine original channel index */
	  ChanIndex = j+RunMode->MaxChans2Add[i]+1;  

	  /* Assign Output Profile to each individual channel */
	  for(k=0;k<RunMode->NBins;k++){
	    TempProfs[RunMode->NumEffChans].rstds[k] = 
	      InputProfs[ChanIndex].rstds[k];
	    TempProfs[RunMode->NumEffChans].rstdq[k] = 
	      InputProfs[ChanIndex].rstdq[k];
	    TempProfs[RunMode->NumEffChans].rstdu[k] = 
	      InputProfs[ChanIndex].rstdu[k];
	    TempProfs[RunMode->NumEffChans].rstdv[k] = 
	      InputProfs[ChanIndex].rstdv[k];
	  }
	  /* Assign central frequencies */
	  OutputFreq[RunMode->NumEffChans] = hdr->obs.ChanFreq[ChanIndex];
	  OutputHead[RunMode->NumEffChans] = HeadLine[ChanIndex];
	  SubOutHdr->DumpRefPhase[RunMode->NumEffChans] = 
	    SubDumpHdr->DumpRefPhase[ChanIndex];
	  SubOutHdr->DumpRefPeriod[RunMode->NumEffChans] =
	    SubDumpHdr->DumpRefPeriod[ChanIndex];
	  RunMode->NumEffChans++;  /* Increase number of effective channels */
	  
	  
	}
      }
    }

    /* Pick out channels not combined: */
    notpicked = (hdr->obs.NChan-1)  - 
      RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1];

    for(j=0;j<notpicked;j++){
      if(RunMode->Verbose){
	printf("Channel %d  OR  %f MHz --> going solo\n",
	       j+RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1]+1,
	       hdr->obs.ChanFreq[j+RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1]+1]);
	fflush(stdout);
      }
      /* Determine original channel index */
      ChanIndex = j+RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1]+1;  
      /* Assign Output Profile to each individual channel */
      for(k=0;k<RunMode->NBins;k++){
	TempProfs[RunMode->NumEffChans].rstds[k]=InputProfs[ChanIndex].rstds[k];
	TempProfs[RunMode->NumEffChans].rstdq[k]=InputProfs[ChanIndex].rstdq[k];
	TempProfs[RunMode->NumEffChans].rstdu[k]=InputProfs[ChanIndex].rstdu[k];
	TempProfs[RunMode->NumEffChans].rstdv[k]=InputProfs[ChanIndex].rstdv[k];
      }
      /* Assign central frequencies */
      OutputFreq[RunMode->NumEffChans] = hdr->obs.ChanFreq[ChanIndex];
      OutputHead[RunMode->NumEffChans] = HeadLine[ChanIndex];
      SubOutHdr->DumpRefPhase[RunMode->NumEffChans] = 
	SubDumpHdr->DumpRefPhase[ChanIndex];
      SubOutHdr->DumpRefPeriod[RunMode->NumEffChans] =
	SubDumpHdr->DumpRefPeriod[ChanIndex];
     
      RunMode->NumEffChans++;  /* Increase number of effective channels */
    }

 
  }
  else{
    /* Just duplicate original input profiles and frequencies to TempProfs */

    RunMode->NumEffChans = hdr->obs.NChan;
    for(i=0;i<RunMode->NumEffChans;i++){
      OutputFreq[i] = hdr->obs.ChanFreq[i];
      OutputHead[i] = HeadLine[i];
      SubOutHdr->DumpRefPhase[i] = SubDumpHdr->DumpRefPhase[i];
      SubOutHdr->DumpRefPeriod[i] = SubDumpHdr->DumpRefPeriod[i];


      for(k=0;k<RunMode->NBins;k++){
	TempProfs[i].rstds[k] = InputProfs[i].rstds[k];
	TempProfs[i].rstdq[k] = InputProfs[i].rstdq[k];
	TempProfs[i].rstdu[k] = InputProfs[i].rstdu[k];
	TempProfs[i].rstdv[k] = InputProfs[i].rstdv[k];
      }
    }
  }

  if(RunMode->Verbose)
    printf("TOTAL NUMBER OF EFFECTIVE CHANNELS = %d\n\n",
	   RunMode->NumEffChans);fflush(stdout);

  /* HERE: Apply bindown (if requested), baseline subtract, and calculation of */
  /*       lin. pol., phi, and phierr    */
  for(i=0;i<RunMode->NumEffChans;i++){
    if(RunMode->BinDown) {
      BinDown(RunMode, &TempProfs[i], &OutputProfs[i]);
    }
    else {
      memcpy(&OutputProfs[i], &TempProfs[i], sizeof(struct StdProfs));
    }
    
    /* If we supply the standard profile then we can correct PA */
/*     if(!RunMode->MakeStd)  */
/*       GetAngle(RunMode, &OutputProfs[i], StdProfile);  */

    MakePol(RunMode, RunMode->NBinsOut, &OutputProfs[i]);


    
  }

  /* Test that output channels are the same as I thought from set up */
  if (test_file){
    Ftest=fopen("test_chans2.dat","a");
    fprintf(Ftest,"Original Channels: \n\n");
    for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
      fprintf(Ftest,"  %6.1lf  ",hdr->obs.ChanFreq[i_chan]);
    }
    fprintf(Ftest,"\n\n\n");
    fprintf(Ftest,"Combined Channels (Total %d):\n\n", RunMode->NumEffChans);
    for(i_chan=0; i_chan<RunMode->NumEffChans; i_chan++){
      fprintf(Ftest,"Out Channel %d: \n", i_chan);
      fprintf(Ftest," %d      --->  %d\n",
	      RunMode->FirstChanAdd[i_chan],RunMode->LastChanAdd[i_chan]);
      fprintf(Ftest," %6.1lf  --->  %6.1lf\n\n", 
	     hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan]], 
	     hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan]]);
      fprintf(Ftest,"Total Scans: %d\n",RunMode->TotScans[nscan*RunMode->NumEffChans + i_chan]);
      fprintf(Ftest,"Total Omit: %d\n\n",RunMode->TotOmit[nscan*RunMode->NumEffChans + i_chan]);
    }
    fclose(Ftest);
  }

  free(AddProfs);
  free(TempProfs);

}
