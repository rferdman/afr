/* Program to organize command line options */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "ASPFitsReader.h"

int GetOptions(struct RunVars *RunMode, struct CalVars *CalMode, 
		Cmdline *Cmd, struct ASPHdr *hdr)
{

  int    i, LastSlashIndex;
  char   TempChar[20];
  int    RootPosition;


  /* Default values */

  /**** Flag Options ****/

  RunMode->Header       = Cmd->HeaderP;
  RunMode->FlipPA       = Cmd->FlipPAP;
  RunMode->MakeRaw      = Cmd->MakeRawP;
  RunMode->Scale        = Cmd->ScaleP;
  RunMode->Verbose      = Cmd->VerboseP;
  RunMode->Cal          = Cmd->CalfileP;
  RunMode->ThetaBBFlag  = Cmd->ThetaBBfileP;
  RunMode->Swap         = Cmd->SwapP;
  RunMode->OldFits      = Cmd->OldFitsP;
  RunMode->NoBase       = Cmd->NoBaseP;
  

  /**** String Options ****/  

  /* Find position of last "." */
  
  RootPosition = -1;
  for (i=strlen(RunMode->Infile)-1;i>0;i--){
    if(!strncmp(&RunMode->Infile[i],".",1)){
      RootPosition = i;
      break;
    }
  }
  if (RootPosition == -1) RootPosition = strlen(RunMode->Infile) - 1;
   
    
/* Mandatory option */
/*  strcpy(RunMode->Source,Cmd->Source); 
  if(strcmp(RunMode->Source, hdr->target.PSRName)!=0){
    printf("Source name does not match that in FITS Header! Exiting...\n");
    fflush(stdout);
    return -1;
  } */

  strcpy(RunMode->Source,hdr->target.PSRName);

  /* if we are just looking at the cal scan file here, make sure that the 
   * Source Name has a "c" in front to distinguish */
  if (strcmp(hdr->gen.ObsMode,"CAL") == 0){
    sprintf(RunMode->Source,"c%s",hdr->target.PSRName);    
  }

  /* If no outfile root is given on the command line, use PSR name and some of  
   * the input file name to create one */
  LastSlashIndex = -1;
  if(!(Cmd->OutfileRootP)) {
    for(i=0;i<strlen(RunMode->Infile);i++){
      if(!strncmp(&RunMode->Infile[i],"/",1))
	LastSlashIndex = i;
    }
    strncpy(TempChar, &RunMode->Infile[LastSlashIndex+1],12);
    strcpy(&TempChar[12],"\0");
    //    sprintf(RunMode->OutfileRoot,"%s.%s",RunMode->Source,TempChar);
    sprintf(RunMode->OutfileRoot,"%s.%s.%s",RunMode->Source,TempChar,hdr->obs.ObsvtyCode);
    if(RunMode->Verbose)
      printf("Infile is %s, OutFileRoot is %s\n",RunMode->Infile, 
	     RunMode->OutfileRoot);fflush(stdout); 
  }
  else {
    strcpy(RunMode->OutfileRoot,Cmd->OutfileRoot);
    if(RunMode->Verbose)
      printf("OutFileRoot is %s\n",RunMode->OutfileRoot);fflush(stdout); 
  }

  /* If no standard profile is given, then don't create one  */
  if(!(Cmd->StdfileP)) {
    printf("No standard profile given.  Will NOT rotate ");
    printf("input profiles to match.\n");fflush(stdout);
    /*    if (strcmp(hdr->gen.ObsMode,"CAL") == 0){
      strncpy(ShortName,RunMode->Source,5);
      ShortName[5] = '\0';
      sprintf(RunMode->Stdfile,"%4s.std",ShortName);
    }
    else{
      strncpy(ShortName,RunMode->Source,4);
     ShortName[4] = '\0';
      sprintf(RunMode->Stdfile,"%4s.std",ShortName);
      } */
  } 
  else {
    strcpy(RunMode->Stdfile,Cmd->Stdfile);
  }    


  if(RunMode->ThetaBBFlag)
    strcpy(RunMode->ThetaBBfile,Cmd->ThetaBBfile);
  
  //printf("THETABBFILE = %s\n ",RunMode->ThetaBBfile);fflush(stdout);



  /*** CAL OPTIONS ***/

  if(RunMode->Cal){ 
    /* If we are calibrating... */

    strcpy(CalMode->Calfile,"\0");
    strcpy(CalMode->Calfile,Cmd->Calfile);


	  /* if we use "default" Mark4 values -- ARECIBO ONLY */

	  /*if(RunMode->FSky > 425. && RunMode->FSky < 435){
	     if(CalMode->CalStrength == 0){
	      CalMode->JyPerCal[0] = 1.041;
	      CalMode->JyPerCal[1] = 1.091;
	    }
	    else if (CalMode->CalStrength == 1){
	      CalMode->JyPerCal[0] = 9.250;    
	      CalMode->JyPerCal[1] = 9.283;
	      } 
	    else{
	      printf("Cal Mode = 0 (low) or 1 (high) only!!!! Exiting...\n");
	      fflush(stdout);
	      return -4;
	    }
	    } */
	  /* else if (RunMode->FSky > 1399. && RunMode->FSky < 1411. 
		   && CalMode->CalStrength == 0) {
	    CalMode->JyPerCal[0] =  0.146;
	    CalMode->JyPerCal[1] =  0.079;
	    } */

  }


  /*** END CAL OPTIONS ***/

  


  /**** Float Options ****/

  //  RunMode->FSky = Cmd->FSky;  /* Mandatory option */
  RunMode->FSky = hdr->obs.FSkyCent;  




  /**** Integer Options ****/
  RunMode->NBins = hdr->redn.RNBinTimeDump;
  RunMode->BinDown = 0; 

  if(!Cmd->NBinsOutP) {
    RunMode->NBinsOut = RunMode->NBins;
  }
  else{
    RunMode->NBinsOut = Cmd->NBinsOut;  /* Output number of bins */
  }
  if(RunMode->NBinsOut < RunMode->NBins) RunMode->BinDown = 1;


  /* Add Dumps together */
  if(Cmd->AddDumpsP && (Cmd->AddDumps > 1)){
    RunMode->AddDumps = Cmd->AddDumps;
    RunMode->NumEffDumps = (int)(ceil(((double)RunMode->NDumps) 
				      / ((double)RunMode->AddDumps)));
    RunMode->NOutDumps = (int)(ceil(((double)RunMode->NDumps) 
				      / ((double)RunMode->AddDumps)));
    if(RunMode->AddDumps > RunMode->NDumps){
      RunMode->AddDumps = RunMode->NDumps;
      printf("You asked to add more dumps together than exist!\n");
      printf("Setting number of dumps to add equal to total number of dumps.\n");
      fflush(stdout);
    }
  } 
  /* option, but no arguments?  why, add all dumps together...! */
  else if(Cmd->AddDumpsP && Cmd->AddDumpsC == 0){ 
    RunMode->AddDumps = RunMode->NDumps;
    RunMode->NumEffDumps = 1;
    RunMode->NOutDumps = 1;
  }
  else{
    /* Default value will be 1, obviously -- no dumps added together */
    RunMode->AddDumps = 1;
    RunMode->NumEffDumps = RunMode->NDumps;
    RunMode->NOutDumps = RunMode->NDumps;
  }

  
#if 0
  /* If user uses a file to list the omissions, then make sure that they 
   * don't use the -scanomit option, and vice-versa */
  if (Cmd->OmitfileP && Cmd->ScanOmitP) {
    printf("Either supply Scan Omissions file or list omissions on command \n");
    printf("line.  Cannot do both! Exiting...\n");fflush(stdout);
    return -10;
  }
#endif

  /* Figure out sideband */
  if (hdr->obs.ChanFreq[0] < hdr->obs.ChanFreq[hdr->obs.NChan-1])
    RunMode->Sideband = 1;   /* upper sideband */
  else
    RunMode->Sideband = -1;  /* lower sideband */
  

  //  printf("MinChan = %lf,  MaxChan = %lf\n\n",MinChan,MaxChan);fflush(stdout);

#if 0
  /* if user has given omissions on command line... 
   * format is "dump freq dump freq ...etc. ..."  */
  if(Cmd->ScanOmitP){
    /* must be an even number of arguments */
    if(Cmd->ScanOmitC%2 != 0){
      printf("-scanomit option must have even number of arguments, with \n");
      printf("format -scanomit dump1 freq1 dump2 freq2 ...\n");fflush(stdout);
      return -12;
    }
    RunMode->NScanOmit = Cmd->ScanOmitC/2;
    /* Now, parse option */
    for (i=0;i<RunMode->NScanOmit;i++){
      RunMode->DumpOmit[i] = Cmd->ScanOmit[2*i];
      FreqOmit[i] = (double)(Cmd->ScanOmit[2*i+1]);
    }
  }
#endif

  /*  if(GetChans(hdr, Cmd, RunMode) < 0) {
    printf("Could not parse channel options.  Exiting...\n");
    return -10;
    } */

  /* Test command line */

  if(RunMode->Verbose) {
    printf("\nCOMMAND LINE OPTIONS:\n");fflush(stdout);
    printf("=====================\n");fflush(stdout);
    printf("Flag options:\n");fflush(stdout);
    printf("-------------\n");fflush(stdout);
    printf("RunMode->Header      = %d\n",RunMode->Header);fflush(stdout);
    printf("RunMode->FlipPA      = %d\n",RunMode->FlipPA);fflush(stdout);
    printf("RunMode->MakeRaw     = %d\n",RunMode->MakeRaw);fflush(stdout);
    printf("RunMode->Scale       = %d\n",RunMode->Scale);fflush(stdout);
    printf("RunMode->Verbose     = %d\n\n",RunMode->Verbose);fflush(stdout);
    printf("String options:\n");fflush(stdout);
    printf("---------------\n");fflush(stdout);
    printf("RunMode->Infile      = %s\n",RunMode->Infile);fflush(stdout);
    printf("RunMode->OutfileRootP = %d\n",Cmd->OutfileRootP);fflush(stdout);
    printf("RunMode->OutfileRoot = %s\n",RunMode->OutfileRoot);fflush(stdout);
    printf("RunMode->Stdfile     = %s\n",RunMode->Stdfile);fflush(stdout);
    printf("RunMode->Source      = %s\n",RunMode->Source);fflush(stdout);
    printf("CalMode->Calfile     = %s\n",CalMode->Calfile);fflush(stdout);
    printf("RunMode->ThetaBBfile = %s\n\n",RunMode->ThetaBBfile);fflush(stdout);
    printf("Float options:\n");fflush(stdout);
    printf("--------------\n");fflush(stdout);
    printf("RunMode->FSky        = %f\n",RunMode->FSky);fflush(stdout);
    printf("RunMode->Cal         = %d\n",RunMode->Cal);fflush(stdout);
    printf("CalMode->CalStrength = %d\n\n",
	   CalMode->CalStrength);fflush(stdout);
    printf("Int options:\n");fflush(stdout);
    printf("------------\n");fflush(stdout);
    printf("RunMode->BinDown     = %d\n",RunMode->BinDown);fflush(stdout);
    printf("RunMode->NBinsOut    = %d\n",RunMode->NBinsOut);fflush(stdout);
    printf("RunMode->AddDumps    = %d\n",RunMode->AddDumps);fflush(stdout);
    printf("RunMode->NumEffDumps = %d\n",RunMode->NumEffDumps);fflush(stdout);
    printf("RunMode->AddChans    = %d\n\n",Cmd->AddChansP);fflush(stdout);
  }



  return 0;
}


