#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPDefs.h"
#include "ASPCommon.h"
#include "ASPFitsReader.h"

int GetOmit(struct ASPHdr *hdr, Cmdline *Cmd, struct RunVars *RunMode)
{


  int    i,j, nomit[NCHMAX];
  int    *NOmitPerScan, nscanomit;
  int    ThisDump, ThisChan;
  //  int    *OmitFlag;
  double FreqOmit[MAXOMIT];
  double MinChan, MaxChan;
  char   OmitLine[128];
  FILE   *Fomit;

  int    i_chan, i_dump, i_omit, i_chan_add, i_dump_add;
  int    OutChan;
  int    MinAddDump, MaxAddDump;
  int    print_to_file=0;
  FILE   *Ftest;

  int DumpOmitSum; //, ChanOmitSum[NCHMAX];
  int DumpOmitFlag[MAXDUMPS];
  int ChanOmitFlag[NCHMAX];

  /* START OF PARSING OF WHICH FREQ/DUMP COMBOS TO OMIT */
  IZero(nomit, NCHMAX);
  RunMode->NScanOmit = 0;
  
  /* make OmitFlag array to denote (Dump,Freq) combinations that are to be omitted */
  //OmitFlag = (int *)malloc(RunMode->NDumps*hdr->obs.NChan*sizeof(int));

  /* Zero it out */
  IZero(&RunMode->OmitFlag[0],RunMode->NDumps*hdr->obs.NChan);
  
  /* Figure out max and min frequencies */
  if (hdr->obs.ChanFreq[0] < hdr->obs.ChanFreq[hdr->obs.NChan-1]){
    MinChan = hdr->obs.ChanFreq[0];
    MaxChan = hdr->obs.ChanFreq[hdr->obs.NChan-1];
  }
  else{
    MinChan = hdr->obs.ChanFreq[hdr->obs.NChan-1];
    MaxChan = hdr->obs.ChanFreq[0];
  }

  /* If user has given file... format is "dump#  freq(MHz) */
  if (Cmd->ZapfileP) {
    if((Fomit = fopen(Cmd->Zapfile, "r")) == NULL){
      printf("Cannot open file %s for reading scan omissions.  Exiting...\n",
	     Cmd->Zapfile);fflush(stdout);
      return -11;
    }
    printf("\nOmission file name:  %s\n",Cmd->Zapfile);
    RunMode->NScanOmit = 0;

    /* Now read the file and construct OmitFlag array */
    while (fgets(OmitLine, 100, Fomit) != NULL){
      sscanf(OmitLine, "%d %lf", &RunMode->DumpOmit[RunMode->NScanOmit],
	     &FreqOmit[RunMode->NScanOmit]);
      /* Give dump number this simple variable name */
      ThisDump = RunMode->DumpOmit[RunMode->NScanOmit];
      /* Chek first whether frequency is negative.  If so, all channel of the 
	 given dump number are flagged for omission */
      if(FreqOmit[RunMode->NScanOmit] < 0.0){
	for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++)
	       RunMode->OmitFlag[ThisDump*hdr->obs.NChan + i_chan] = 1;
      }
      else {
      /* Now we know frequency is not negative. So, map frequency to integer 
	 channel array number */
	if((ThisChan = Freq2Chan(FreqOmit[RunMode->NScanOmit], 
				 hdr->obs.ChanFreq, hdr->obs.NChan)) < 0) {
	  /* If mapping unsuccessful, quit out with error */
	  fprintf(stderr, "Error in scan omission file. Frequency %f does not ", 
		  FreqOmit[RunMode->NScanOmit]);
	  fprintf(stderr, "map to any channel for this data file. \n");
	  return -12;
	}
	
	/* Now check whether dump number is negative. If so, all dumps at the 
	   given frequency are flagged for omission */
	if(ThisDump < 0){
	  for (i_dump=0; i_dump<RunMode->NDumps ; i_dump++)
	    	  RunMode->OmitFlag[i_dump*hdr->obs.NChan + ThisChan] = 1;
	}
	else{
	  /* Finally, if we got to here, then neither frequency nor dump number 
	     are negative.  In this case, flag this dump/channel combination 
	     as one for omission */
	  RunMode->OmitFlag[ThisDump*hdr->obs.NChan + ThisChan] = 1;
	  RunMode->NScanOmit++;
	  /*	  printf("OmitFlag index = %d, NScanOmit = %d\n\n", 
		 ThisDump*hdr->obs.NChan + ThisChan, 
		 RunMode->NScanOmit);fflush(stdout);  */
	}

      }
    }

    fclose(Fomit);
  }


  /**** print out OmitFlag matrix ****/
  if (print_to_file) {
    Ftest=fopen("test_omit.dat", "w");
    fprintf(Ftest, "  ");
    for (i_chan=0; i_chan<hdr->obs.NChan; i_chan++)
      fprintf(Ftest, " %6.1lf ",hdr->obs.ChanFreq[i_chan]);
    fprintf(Ftest, "\n\n");
 
    for (i_dump=0; i_dump<RunMode->NDumps; i_dump++){
      fprintf(Ftest, "%d  ",i_dump);
      for (i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
	fprintf(Ftest,"   %d    ",RunMode->OmitFlag[i_dump*hdr->obs.NChan + i_chan]);
      }
      fprintf(Ftest, "\n");
    }
    fclose(Ftest);
  }

  /////  Maybe can get rid of this soon /////
  if(Cmd->ScanOmitP){
    RunMode->NumAllDumpOmit = Cmd->ScanOmitC;
    for (i=0;i<RunMode->NumAllDumpOmit;i++) { 
      RunMode->AllDumpOmit[i]=Cmd->ScanOmit[i];
    }   
  }
  else {
    RunMode->NumAllDumpOmit = 0;
  }
  ///////////////////////////////////////////

  if(RunMode->Verbose){
    printf("\nOMISSIONS FROM INPUT FILE %s:\n\n",RunMode->Infile);
    printf("Dump   Frequency\n");fflush(stdout);
  }
    /* Check to make sure that values given in file are within current
     * scan's freqs and dumps, then see if all frequencies within a 
     * scan are chosen, in which case add that scan to the AllDumpOmit
     * array */

  /* First, create array to store channel omission info for each scan */

  NOmitPerScan=(int *)malloc(RunMode->NDumps*hdr->obs.NChan*sizeof(int));
  IZero(NOmitPerScan, RunMode->NDumps*hdr->obs.NChan);

  for(i=0;i<RunMode->NScanOmit;i++){
    if(RunMode->Verbose)
      printf("%4d   %9lf\n",
	     RunMode->DumpOmit[i],FreqOmit[i]);
    
    if(RunMode->DumpOmit[i] > (RunMode->NDumps-1)){
      printf("Individual dump numbers must be less than\n");
      printf("total number of dumps (= %d for the current file) \n",
	     RunMode->NDumps);fflush(stdout);
      return -13;
    }
    if(FreqOmit[i] < MinChan || FreqOmit[i] > MaxChan){
      printf("Frequency omission out of range (%lf -> %lf)\n",
	     MinChan, MaxChan);fflush(stdout);
      return -14;
    }

    ///////  can probably make this obsolete soon ////////

    /* Figure out which channel corresponds to each frequency */
    RunMode->ChanOmit[i] = -1;
    for (j=0;j<hdr->obs.NChan;j++){
      if (fabs(FreqOmit[i] - hdr->obs.ChanFreq[j]) <= DBLEPS){
	//      if (FreqOmit[i] == hdr->obs.ChanFreq[j]){
	RunMode->ChanOmit[i] = j;
	nomit[j]++;
	/* Include this channel as an omission for this dump */
        NOmitPerScan[RunMode->DumpOmit[i]*hdr->obs.NChan + j] = 1;
	break;
      }	
    }
    if(RunMode->ChanOmit[i] < 0){
      printf("Frequency %lf is not found in data. Exiting...\n",
	     FreqOmit[i]);fflush(stdout);
      return -15;
    }
    ///////////////////////////////////////////////////////


  }
  

  /////////  can probably make this obsolete soon  ////////

  /* Figure out which dumps have all frequencies omitted */
  for(i=0;i<RunMode->NDumps;i++){
    nscanomit=0;
    for(j=0;j<hdr->obs.NChan;j++){
      nscanomit += NOmitPerScan[i*hdr->obs.NChan + j];
    }
    if (nscanomit == hdr->obs.NChan) {
      RunMode->AllDumpOmit[RunMode->NumAllDumpOmit++]=i;
      /*printf("Dump %d added to list to be omitted entirely\n\n"
	,i);fflush(stdout); */
    }
  }
  /////////////////////////////////////////////////////////


  /* do same as above but using the OmitFlag matrix */
  IZero(DumpOmitFlag, MAXDUMPS);
  
  //  IZero(ChanOmitSum, NCHMAX);

  /* Firstly, identify scan #'s for which all chans are flagged to be omitted */
  for(i_dump=0; i_dump<RunMode->NDumps; i_dump++) {
    DumpOmitSum = 0;
    for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
      DumpOmitSum += RunMode->OmitFlag[i_dump*hdr->obs.NChan + i_chan];    
    }
    if(DumpOmitSum == hdr->obs.NChan) {
      printf("Dump %d added to list to be omitted entirely.\n\n", i_dump);
      DumpOmitFlag[i_dump]=1;
    }
  }
  
  /* Now, flag those scans for omission given on command line */
  for (i_omit=0; i_omit<Cmd->ScanOmitC; i_omit++){
    DumpOmitFlag[ Cmd->ScanOmit[i_omit] ] = 1;
    /* Copy 1's to all channels in the OmitFlag array */
    for (i_chan=0;i_chan<hdr->obs.NChan; i_chan++)
      RunMode->OmitFlag[Cmd->ScanOmit[i_omit]*hdr->obs.NChan + i_chan]=1;
  }

  /* NOW, all scans to be omitted entirely are recorded in OmitFlag AND 
     DumpOmitFlag  */


  /* print DumpOmitFlag: */
  /*  for(i_dump=0; i_dump<RunMode->NDumps; i_dump++) 
      printf("%d    %d\n",i_dump, DumpOmitFlag[i_dump]); */
  
  
  
  /* END OF PARSING OF WHICH FREQ/DUMP COMBOS TO OMIT */
  
  /* START OF PARSING OF WHICH FREQ CHANNEL TO OMIT COMPLETELY */

  if(Cmd->FreqOmitP || Cmd->ChanOmitP){
    printf("All scans at the following frequencies will be omitted:\n");
    
    //// maybe gone soon ////
    for(j=0;j<hdr->obs.NChan;j++)
      RunMode->ZapChan[j]=0;
    ////////////////////////
  }
  
  if (Cmd->FreqOmitP){
    
    IZero(ChanOmitFlag,NCHMAX);
    
    for(i=0;i<Cmd->FreqOmitC;i++){
      
      if(Cmd->FreqOmit[i] < MinChan || Cmd->FreqOmit[i]> MaxChan){
	printf("Frequency omission (%lf) out of range (%lf -> %lf)\n",
	       Cmd->FreqOmit[i], MinChan, MaxChan);fflush(stdout);
	return -16;
      }
      
      /* Map frequency to channel number */
      if((ThisChan = Freq2Chan(Cmd->FreqOmit[i], 
			       hdr->obs.ChanFreq, hdr->obs.NChan)) < 0) {
	printf("Error in scan omission file.\n");
	return -12;
      }
      /* Flag it */
      ChanOmitFlag[ThisChan]=1;
      /* Now add flags to OmitFlag matrix for all dumps at this channel */
      for (i_dump=0; i_dump<RunMode->NDumps; i_dump++)
	RunMode->OmitFlag[i_dump*hdr->obs.NChan + ThisChan]=1;

      
      //// maybe gone soon ////
      
      /* Flag certain channels if they have been chosen for omission */
      for(j=0;j<hdr->obs.NChan;j++){
	if(hdr->obs.ChanFreq[j] == Cmd->FreqOmit[i]){
	  RunMode->ZapChan[j]=1;
	  break;
	}
      }
      /////////////////////////
   
      printf("Channel %4d = %.3lf MHz\n", 
	     ThisChan, hdr->obs.ChanFreq[ThisChan]);
         
    }


    if (print_to_file) {
      Ftest=fopen("test_omit2.dat", "w");
      fprintf(Ftest, "Channels to be omitted: \n\n");
      for(i=0;i<Cmd->FreqOmitC;i++)
	fprintf(Ftest, "%6.1lf   ", Cmd->FreqOmit[i]);
      fprintf(Ftest, "\n\n\n");

      for (i_chan=0; i_chan<hdr->obs.NChan; i_chan++)
	fprintf(Ftest, " %6.1lf ",hdr->obs.ChanFreq[i_chan]);
      fprintf(Ftest, "\n\n");
      
      for (i_dump=0; i_dump<RunMode->NDumps; i_dump++){
	fprintf(Ftest, "%4d  ",i_dump);
	for (i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
	  fprintf(Ftest,"   %d    ",RunMode->OmitFlag[i_dump*hdr->obs.NChan + i_chan]);
	}
	fprintf(Ftest, "\n");
      }
      fclose(Ftest);
    }
    
    
    /* Also, if we have manually omitted all scans of a given frequency above,
       then just omit frequency */
    
    //// maybe gone soon ////
    for(j=0;j<hdr->obs.NChan;j++){
      if (nomit[j] >= RunMode->NDumps)
	RunMode->ZapChan[j]=1;
    }
    /////////////////////////
    
 
  }


  /* Now do the same for the ChanOmit command line option */

  /* This takes integer channel number and adds all dumps in that channel 
     to omit flag array */

  if (Cmd->ChanOmitP){

    for(i=0;i<Cmd->ChanOmitC;i++){

      if(abs(Cmd->ChanOmit[i]) > hdr->obs.NChan -1){
	printf("Frequency channel omission (%d) out of range ",
	       Cmd->ChanOmit[i]);
	printf("(channel number 0 -> %d)\n", hdr->obs.NChan-1);fflush(stdout);
	return -16;
      }
      
      
      if (Cmd->ChanOmit[i] < 0) 
	ThisChan = hdr->obs.NChan + Cmd->ChanOmit[i];
      else
	ThisChan = Cmd->ChanOmit[i];

      for (i_dump=0; i_dump<RunMode->NDumps; i_dump++){
	RunMode->OmitFlag[i_dump*hdr->obs.NChan + ThisChan] = 1;
	////// Following is pretty much obsolete //////
	RunMode->ZapChan[ThisChan] = 1;
      }
      
      printf("Channel %4d = %.3lf MHz\n", 
	     ThisChan, hdr->obs.ChanFreq[ThisChan]);

    }

  }

  printf("\n");
  
  /* END OF PARSING OF WHICH FREQ CHANNEL TO OMIT COMPLETELY */


  /* Figure out output dump and channel structure and calculate denominators
     for averaging by taking into account omitted scans */

  /* Initialize */
  IZero(RunMode->TotScans, MAXDUMPS*NCHMAX);
  IZero(RunMode->TotOmit, MAXDUMPS*NCHMAX);

  /* Do one (output) channel at a time */
  for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
    /* Add OmitFlags in dumps, and increment by number of consecutive 
       dumps to add together */
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
      
      /* Now loop over beginning and end channels of each output channel */
      /* Note that we put "<=" here since we want the range to be inclusive */
      for (i_chan_add=RunMode->FirstChanAdd[i_chan]; 
	   i_chan_add<=RunMode->LastChanAdd[i_chan]; i_chan_add++){ 
	
      /* Now loop over beginning and end dumps of each output dump */
	for (i_dump_add=MinAddDump; i_dump_add<MaxAddDump; i_dump_add++){
	  
	  // printf("i_chan_add= %d,   i_dump_add = %d\n", i_chan_add,i_dump_add);
	  RunMode->TotScans[i_dump*RunMode->NOutChans + i_chan]++;
	  RunMode->TotOmit[i_dump*RunMode->NOutChans + i_chan] += 
	    RunMode->OmitFlag[i_dump_add*hdr->obs.NChan + i_chan_add];
	  
	}
	
      }
      
    }

  }
  
  if (print_to_file) {
    Ftest=fopen("test_omit_totals.dat", "w");
    fprintf(Ftest, "Original NDumps = %d;  AddDumps = %d\n\n",
	    RunMode->NDumps, RunMode->AddDumps);
    fprintf(Ftest, "  Total Input Scans in each Output Dump/Channel combination:\n\n");
    fprintf(Ftest, "   ");
   

    for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
      /* OutChan is average channel */
      OutChan = (int)(((double)(RunMode->FirstChanAdd[i_chan] + 
				RunMode->LastChanAdd[i_chan]))/2.);
      fprintf(Ftest, " %6.1lf ",hdr->obs.ChanFreq[OutChan]);
    }
    fprintf(Ftest, "\n\n");
    
    for (i_dump=0; i_dump<RunMode->NOutDumps; i_dump++){
      fprintf(Ftest, "%d  ",i_dump);
      for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
	fprintf(Ftest,"   %d    ",
		RunMode->TotScans[i_dump*RunMode->NOutChans + i_chan]);
      }
      fprintf(Ftest, "\n");
    }
    fprintf(Ftest, "\n\n\n");
    fprintf(Ftest, "  Total Omissions from each Output Dump/Channel combination:\n\n");
    for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
      /* OutChan is average channel */
      OutChan = (int)(((double)(RunMode->FirstChanAdd[i_chan] + 
				RunMode->LastChanAdd[i_chan]))/2.);
      fprintf(Ftest, " %6.1lf ",hdr->obs.ChanFreq[OutChan]);
    }
    fprintf(Ftest, "\n\n");
    
    for (i_dump=0; i_dump<RunMode->NOutDumps; i_dump++){
      fprintf(Ftest, "%d  ",i_dump);
      for (i_chan=0; i_chan<RunMode->NOutChans; i_chan++){
	fprintf(Ftest,"   %d    ",
		RunMode->TotOmit[i_dump*RunMode->NOutChans + i_chan]);
      }
      fprintf(Ftest, "\n");
    }



    fclose(Ftest);
  }
  
  //  free(OmitFlag);


  



   return 0; 


}



