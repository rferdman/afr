#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "ASPFitsReader.h"

/*  Get, parse, and organize channel output options */
int GetChans(struct ASPHdr *hdr, Cmdline *Cmd, struct RunVars *RunMode)
{

  int    i, j, i_chan, i_chan_in, i_chan_out, ChanIndex;
  int    chan_up=0, chan_down=0, chan_hi=0, chan_lo=0, chan_start=-1;
  int    notpicked; 
  int    test_chan_file=0;
  int    NChanSub;
  double BW;
  FILE *Ftest;

  int testarray[10];

  BW = DSum(hdr->obs.ChanWidth, hdr->obs.NChan);

  /* START OF CHANNEL ADDING OPTION PARSING */
  /* Check to see if there are chans to add: */
  RunMode->AddChans = Cmd->AddChansP;
  


  if(Cmd->AddChansP) {

    if(Cmd->NSubsP || Cmd->SubBW || Cmd->SubRefP) {
      printf("Cannot use -nsubs, -subbw, or -subref when using -addchans. ");
      printf("Exiting...\n");
      return -19;
    }

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
  /* Subband options */
  else if(Cmd->NSubsP || Cmd->SubBWP || Cmd->NChanSubP || Cmd->SubRefP) {
    /* Make sure that the user did not combine subband channel-finding
       options */
    if ( (Cmd->NSubsP && Cmd->SubBWP) ||
	 (Cmd->NSubsP && Cmd->NChanSubP) ||
	 (Cmd->SubBWP && Cmd->NChanSubP) ){
      printf("Can only use one of -nchansub, -subbw, or -nsubs at a time. ");
      printf("Exiting...\n");
      return -19;
    }
    /* Make sure that NSubs is > 0 and < NChans */
    if(Cmd->NSubsP && (Cmd->NSubs < 1) ) {
      printf("Argument to -nsubs must be greater than zero.  Exiting...\n");
      return -19;
    }
    if(Cmd->NSubsP && (Cmd->NSubs > hdr->obs.NChan) ) {
      printf("Argument to -nsubs cannot be more than number of channels in \n");
      printf("data file (%d).  Exiting...\n",hdr->obs.NChan);
      return -19;
    }
    /* Make sure that NChanSub is > 0 and < NChans */
    if(Cmd->NChanSubP && (Cmd->NChanSub < 1) ) {
      printf("Argument to -nchansub must be greater than zero.  Exiting...\n");
      return -19;
    }
    if(Cmd->NChanSubP && (Cmd->NChanSub > hdr->obs.NChan) ) {
      printf("Argument to -nchansub cannot be more than number of channels \n");
      printf("in data file (%d).  Exiting...\n",hdr->obs.NChan);
      return -19;
    }
    /* Make sure that if -subref is used, then one of -nsubs -subbw, 
       or -nchansub is also used */
    if(Cmd->SubRefP && !(Cmd->NSubsP || Cmd->SubBWP || Cmd->NChanSubP)){
      printf("Must use one of the -nsubs, -subbw, or -nchansub options when ");
      printf("using the -subref option to fix a reference subband ");
      printf("frequency. Exiting...\n");
      fflush(stdout);
      return -19;
    }

    /* Now calculate the individual subband bandwidth, 
       in number of channels (i.e. not in MHz).  This is done by rounding off 
       (total channels)/(nsubs) to nearest integer number of channels */
    printf("%d input channels\n", hdr->obs.NChan);
    if(Cmd->NSubsP) {
      NChanSub = (int)(floor((double)hdr->obs.NChan/(double)Cmd->NSubs));
      RunMode->NOutChans = Cmd->NSubs;
      printf("%d subbands chosen.\n", Cmd->NSubs);
      printf("%.1lf MHz subband width\n", NChanSub*(BW/(double)hdr->obs.NChan));
    }
    else if(Cmd->SubBWP){ 
      /* Here, round off to nearest multiple of actual input channel BW */
      NChanSub = (int)(floor(Cmd->SubBW/(BW/(double)hdr->obs.NChan) +0.5));
      RunMode->NOutChans = (int)(ceil((double)hdr->obs.NChan/(double)NChanSub));
      printf("%.1lf-MHz subbands chosen.\n", Cmd->SubBW);
    }
    else { // Cmd->NChanSub:
      NChanSub = Cmd->NChanSub;
      RunMode->NOutChans = (int)(ceil((double)hdr->obs.NChan/(double)NChanSub));
      printf("%.1lf MHz subband width\n", NChanSub*(BW/(double)hdr->obs.NChan));
    }
    printf("%d channels per subband.\n",NChanSub);


    /* If user has chosen a reference frequency for subbanding: */
    if(Cmd->SubRefP) {
      /* Make sure that -subref argument corresponds to the boundary of 
	 two channels */
      chan_hi=0;
      chan_lo=0;
      for(i_chan_in=0; i_chan_in<hdr->obs.NChan; i_chan_in++){
	if (Cmd->SubRef == hdr->obs.ChanFreq[i_chan_in] - 2.)
	  chan_hi++;
	if (Cmd->SubRef == hdr->obs.ChanFreq[i_chan_in] + 2.)
	  chan_lo++;
      }

      chan_up = chan_down = 0;
      if(chan_hi + chan_lo == 0) { // bad choice of ref frequency
	printf("Argument to -subref MUST be a frequency (in MHz) that ");
	printf("borders two existing channel frequencies.  The following are ");
	printf("the centre frequencies of the channels in the current data ");
	printf("file:\n\n");
	for(i_chan_in=0; i_chan_in<hdr->obs.NChan; i_chan_in++){
	  printf("%.1lf   ",hdr->obs.ChanFreq[i_chan_in]);
	}
	printf("\n\nAnd thus these are the choices of -subref arguments:\n\n");
	printf("%.1lf   ",hdr->obs.ChanFreq[0]-
	       RunMode->Sideband*2.);
	for(i_chan_in=0; i_chan_in<hdr->obs.NChan; i_chan_in++){
	  printf("%.1lf   ",hdr->obs.ChanFreq[i_chan_in]+
		 RunMode->Sideband*2.);
	}
	/*	printf("%.1lf   ",hdr->obs.ChanFreq[hdr->obs.NChan-1]+
		RunMode->Sideband*2.); */
	printf("\n\n");
	return -19;
      }
      /* Get start channel */
      else if(chan_hi + chan_lo == 1){  /* We are at one of the end frequecies */
	chan_start = Freq2Chan(Cmd->SubRef + 2., 
			      hdr->obs.ChanFreq, hdr->obs.NChan);
	if (chan_start < 0) { /* oops, went wrong direction */
	  chan_start = Freq2Chan(Cmd->SubRef - 2., 
			      hdr->obs.ChanFreq, hdr->obs.NChan);
	}
	/* Now find out if user chose 0th or NChan-th channel */
	if (chan_start == 0) chan_up=1;
	else /* (chan_start == hdr->obs.NChan-1) */ chan_down=1;
      }
      else{ /* channel lands somewhere in the middle */
	/* Choose channel with  next lowest centre frequency to start */
	chan_start = Freq2Chan(Cmd->SubRef - 2., 
			       hdr->obs.ChanFreq, hdr->obs.NChan);
      }

  /*  printf("chan_up = %d, chan_down = %d, chan_hi = %d, chan_lo = %d\n",
	     chan_up, chan_down, chan_hi, chan_lo);
      printf("SubRef = %.1lf,  chan_start = %d\n", Cmd->SubRef, chan_start); */

      /* Now go through two loops:  One going down in channel index and 
	 one going up in channel number.  Start at chan_start. */
      
      /* Reset NOutChans and count as we go along, just in case there is some
	 issue with how the subbands are distributed */
      RunMode->NOutChans=0;
      /* Going down...*/
      for (i_chan_in=chan_start; i_chan_in>=0; i_chan_in-=NChanSub){
      /* Ensure that first output channel begins at input channel zero */
	if(i_chan_in-NChanSub+1 < 0)
	  RunMode->FirstChanAdd[RunMode->NOutChans]=0;
	else
	  RunMode->FirstChanAdd[RunMode->NOutChans]=i_chan_in-NChanSub+1;
	RunMode->LastChanAdd[RunMode->NOutChans]=i_chan_in;	  
 	RunMode->NOutChans++;
      }
      /* Going up... start at the next channel index up from chan_start 
	 since we started from chan_start previously */
      for (i_chan_in=chan_start+1;i_chan_in<hdr->obs.NChan;i_chan_in+=NChanSub){
	RunMode->FirstChanAdd[RunMode->NOutChans]=i_chan_in;
	/* Ensure that last output channel end at final input channel */
	if(i_chan_in+NChanSub-1 > hdr->obs.NChan-1)
	  RunMode->LastChanAdd[RunMode->NOutChans]=hdr->obs.NChan-1;
	else
	  RunMode->LastChanAdd[RunMode->NOutChans]=i_chan_in+NChanSub-1;
	RunMode->NOutChans++;
      }

      /* Now sort First/LastChanAdd arrays.  A straight sort of each should 
	 do since there should not be any overlap between a LastChanAdd and a
	 FirstChanAdd of another set */

      if (Cmd->VerboseP){
	printf("SUBREF PRE-SORT:\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%d   ",RunMode->FirstChanAdd[i_chan_out]);
	printf("\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%.0lf ",hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan_out]]);
	printf("\n\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%d   ",RunMode->LastChanAdd[i_chan_out]);
	printf("\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%.0lf ",hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan_out]]);
	printf("\n\n");
      }
      
      ISort(RunMode->NOutChans, &RunMode->FirstChanAdd[0]);
      ISort(RunMode->NOutChans, &RunMode->LastChanAdd[0]);
      
      
      if (Cmd->VerboseP) {
	printf("SUBREF POST-SORT:\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%d   ",RunMode->FirstChanAdd[i_chan_out]);
	printf("\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%.0lf ",hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan_out]]);
	printf("\n\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%d   ",RunMode->LastChanAdd[i_chan_out]);
	printf("\n");
	for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
	  printf("%.0lf ",hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan_out]]);
	printf("\n\n");
      }
      
      //      exit(0);
    }
    else {
      /* If no reference frequency chosen, start subbanding from channel zero */
      /* Now fill in FirstChanAdd and LastChanAdd, so that new subbanded data 
	 channels are contiguous */
      for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
	RunMode->FirstChanAdd[i_chan_out] = i_chan_out*NChanSub;
	if( ((i_chan_out+1)*NChanSub - 1) > hdr->obs.NChan-1)
	  RunMode->LastChanAdd[i_chan_out] = hdr->obs.NChan-1;
	else
	  RunMode->LastChanAdd[i_chan_out] = (i_chan_out+1)*NChanSub - 1;
      }
      /* Finally, check that last output channel extends to include last 
	 input channel*/
      RunMode->LastChanAdd[RunMode->NOutChans-1] = hdr->obs.NChan-1;
    }
    

    printf("%d output subbbands over %.1lf MHz bandwidth:\n\n", 
	   RunMode->NOutChans, BW);
    printf("START   ");
    for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
      printf("%6.1lf   ", hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan_out]]);
    printf("\nEND     ");
    for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++)
      printf("%6.1lf   ", hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan_out]]);
    printf("\n\n");

    //  exit(0);
   

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
