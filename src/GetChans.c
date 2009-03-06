#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "ASPFitsReader.h"

/*  Get, parse, and organize channel output options */
int GetChans(struct ASPHdr *hdr, Cmdline *Cmd, struct RunVars *RunMode)
{

  int  i, j, i_chan, ChanIndex;
  int  notpicked; 
  int  test_chan_file=0;
  FILE *Ftest;

  /* START OF CHANNEL ADDING OPTION PARSING */
  /* Check to see if there are chans to add: */
  RunMode->AddChans = Cmd->AddChansP;
   
  if(Cmd->AddChansP) {
    if(Cmd->AddChansC%2 > 0) {  /* if not even number of args */
      printf("-addchans:  must have even number of arguments (currently you have %d)! Exiting...\n",Cmd->AddChansC);
      fflush(stdout);
      return -17;
    }
    
    RunMode->SetsofChans2Add = Cmd->AddChansC/2;  /* self explanatory */

    if(RunMode->Verbose)
      printf("SetsofChans2Add = %d\n",RunMode->SetsofChans2Add);fflush(stdout);

    if(RunMode->SetsofChans2Add == 0) {

      /* If no args, then all channels are added together */
      RunMode->SetsofChans2Add = 1;
      RunMode->MinChans2Add[0] = 0;
      RunMode->MaxChans2Add[0] = hdr->obs.NChan-1;
      RunMode->NumChans2Add[0] = hdr->obs.NChan;
    }
    else {
    
      /* Organize which ranges of channels to add together */
      for(i=0;i<RunMode->SetsofChans2Add;i++) {


	/* Here: set values for Min/MaxChans2Add[i] based on subband 
	   preferences*/

	RunMode->MinChans2Add[i] = Cmd->AddChans[2*i];
	RunMode->MaxChans2Add[i] = Cmd->AddChans[2*i+1];
	RunMode->NumChans2Add[i] = RunMode->MaxChans2Add[i] - 
                                   RunMode->MinChans2Add[i] + 1;

	/* make sure min < max */
	if(RunMode->MinChans2Add[i] >= RunMode->MaxChans2Add[i]) {
	  printf("-addchans: Bottom of range (%d) must be smaller\n",
		 RunMode->MinChans2Add[i]);
	  printf("           than top (%d)! Exiting...\n",
		 RunMode->MaxChans2Add[i]);
	  fflush(stdout);
	  return -18;
	}
	/* make sure min and/or max < number of channels */
	if((RunMode->MinChans2Add[i] > (hdr->obs.NChan-1)) || 
	   (RunMode->MaxChans2Add[i] > (hdr->obs.NChan-1))) {
	  printf("-add: There are only %d channels! Exiting...\n",
		 hdr->obs.NChan);fflush(stdout);
	  return -19;
	}
	/* make sure ranges are in increasing order and not overlapping */
	if(RunMode->SetsofChans2Add > 1 && i > 0){
	  if(RunMode->MinChans2Add[i] <= RunMode->MinChans2Add[i-1] || 
	     RunMode->MaxChans2Add[i] <= RunMode->MaxChans2Add[i-1] ||
	     RunMode->MinChans2Add[i] <= RunMode->MaxChans2Add[i-1] ){
	    printf("-add: Channel ranges must be in increasing order\n");
	    printf("      and not overlapping. Exiting...\n");fflush(stdout);
	    return -20;
	  }

	}
	    
      }  /* END for(i=0;i<RunMode->SetsofChans2Add;i++) */

    }

  }  /*  if(RunMode->AddChans) */
  RunMode->NumEffChans = 0;  /* Number of "effective" channels */
  /* END OF ADDING CHANNELS OPTION PARSING */
  


  /* Find NumEffChans */
  RunMode->NOutChans = 0;

  if(Cmd->AddChansP){
    if(RunMode->Verbose) printf("\nCHANNELS TO ADD:\n");
    /* Pick out channels not combined: */
    notpicked = RunMode->MinChans2Add[0];
    RunMode->NumEffChans += notpicked;

    for(j=0;j<notpicked;j++){
      if(RunMode->Verbose) 
	printf("Channel %d  OR  %f MHz --> going solo\n",
	       j,hdr->obs.ChanFreq[j]);
      

      /* Begin to set up arrays of first and last channels to add together */

      /* for those channels not added together, first == last */
      RunMode->FirstChanAdd[RunMode->NOutChans] = j;
      RunMode->LastChanAdd[RunMode->NOutChans]  = j;
      RunMode->NOutChans++;


    }
    for(i=0;i<RunMode->SetsofChans2Add;i++){
      if(RunMode->Verbose)
	printf("Channels %d to %d  OR  %f to %f MHz\n",
	       RunMode->MinChans2Add[i],RunMode->MaxChans2Add[i],
	       hdr->obs.ChanFreq[RunMode->MinChans2Add[i]],
	       hdr->obs.ChanFreq[RunMode->MaxChans2Add[i]]);
	
      /* Take note of first and last channels for current output channel */
      RunMode->FirstChanAdd[RunMode->NOutChans] = RunMode->MinChans2Add[i];
      RunMode->LastChanAdd[RunMode->NOutChans]  = RunMode->MaxChans2Add[i];
      RunMode->NOutChans++;

      /* Pick out channels not combined: */
      if(i<RunMode->SetsofChans2Add-1){
	notpicked = RunMode->MinChans2Add[i+1] - RunMode->MaxChans2Add[i] - 1;
	RunMode->NumEffChans += notpicked;

  
	for(j=0;j<notpicked;j++){

	  ChanIndex = j+RunMode->MaxChans2Add[i]+1;

	  if(RunMode->Verbose){
	    printf("Channel %d  OR  %f MHz --> going solo\n",
		   ChanIndex, hdr->obs.ChanFreq[ChanIndex]);
	  }
	  
	  /* Take note of output channels for not-added channels */
	  RunMode->FirstChanAdd[RunMode->NOutChans] = ChanIndex;
	  RunMode->LastChanAdd[RunMode->NOutChans]  = ChanIndex;
	  RunMode->NOutChans++;

	}
      }
    }
    RunMode->NumEffChans += RunMode->SetsofChans2Add;

    /* Pick out channels not combined: */
    notpicked = (hdr->obs.NChan-1)  - 
                RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1];
    RunMode->NumEffChans += notpicked;

    for(j=0;j<notpicked;j++){
      
      ChanIndex = j+RunMode->MaxChans2Add[RunMode->SetsofChans2Add-1]+1;

      if(RunMode->Verbose)
        printf("Channel %d  OR  %f MHz --> going solo\n",
	       ChanIndex, hdr->obs.ChanFreq[ChanIndex]);
      
      RunMode->FirstChanAdd[RunMode->NOutChans] = ChanIndex;
      RunMode->LastChanAdd[RunMode->NOutChans]  = ChanIndex;
      RunMode->NOutChans++;
      
    }
  }
  else{
    /* Channels transfer exactly if we are not adding channels */
    for(i_chan=0;i_chan<hdr->obs.NChan;i_chan++){ 
      RunMode->FirstChanAdd[i_chan] = i_chan;
      RunMode->LastChanAdd[i_chan]  = i_chan;
    }
    RunMode->NumEffChans = hdr->obs.NChan;
    RunMode->NOutChans = hdr->obs.NChan;
  }
  if(RunMode->Verbose) {
    printf("TOTAL NUMBER OF EFFECTIVE CHANNELS = %d\n\n",
	   RunMode->NumEffChans);fflush(stdout);
  }


  /*  test file to check chan assignments */
  if (test_chan_file) {
    Ftest=fopen("test_chan_file.dat","w");
    fprintf(Ftest,"Original Channels: \n\n");
    for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
      fprintf(Ftest,"  %6.1lf  ",hdr->obs.ChanFreq[i_chan]);
    }
    fprintf(Ftest,"\n\n\n");
    fprintf(Ftest,"Combined Channels (Total %d):\n\n", RunMode->NOutChans);
    for(i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
      fprintf(Ftest,"Out Channel %d: \n", i_chan);
      fprintf(Ftest," %d      --->  %d\n",
	      RunMode->FirstChanAdd[i_chan],RunMode->LastChanAdd[i_chan]);
      fprintf(Ftest," %6.1lf  --->  %6.1lf\n\n", 
	     hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan]], 
	     hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan]]);
    }
    fclose(Ftest);
  }

  return 1;

}
