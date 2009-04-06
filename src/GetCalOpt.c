/* Program to organize command line options */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "CalCmdLine.h"
#include "ASPCal.h"

int GetCalOpt(struct RunVars *RunMode, struct CalVars *CalMode, 
	      Cmdline *CalCmd, struct ASPHdr *CalHdr)
{

  int i, pol;
  int LastSlashIndex;
  int RootPosition;
  char TempChar[20];

  /**** Flag Options ****/

  RunMode->FlipPA   = CalCmd->FlipPAP;
  RunMode->MakeRaw  = CalCmd->MakeRawP;
  RunMode->Verbose  = CalCmd->VerboseP; 
  RunMode->AddDumps = CalCmd->AddDumpsP;
  RunMode->OldFits = 0;

  /* Make sure that user has given Tcal if Contfile is not given */
  if (!CalCmd->ContfileP) {
    if (!CalCmd->TcalP || !CalCmd->GainP){
      printf("Must give Tcal and telescope Gain if Continuum Cal scan data\n");
      printf("file is not given. Exiting...\n");fflush(stdout);
      return -1;
    }
    if(CalCmd->ChooseMethodP || CalCmd->ConstTsysP || CalCmd->NBadThreshP){
      printf("The options -choose, -tsys, and -nbad are not available if ");
      printf("we are not using continuum calibration mode. They will thus be");
      printf("ignored here...\n");
      fflush(stdout);
      //      return -1;
    }
  }
  else{ // i.e. if we have chosen to use continuum cals
    /* Make sure that cal method flags don't cross */
    if(CalCmd->ConstTsysP && CalCmd->ChooseMethodP) {
      printf("Please choose EITHER constant Tsys method OR having ASPCal ");
      printf("choose method automatically.  Exiting...\n");
      fflush(stdout);
      return -1;
    }
    
    if(CalCmd->ChooseMethodP) {
      printf("ASPCal will choose the method of Cal Flux calculation.\n");
      if(CalCmd->ChooseMethodC == 0){ // i.e. no arguments --> default value
	CalCmd->ChooseMethod = 15.0;  //default value
	printf("Will use the default minimum threshold of %5.1f%% difference ",
	       CalCmd->ChooseMethod);
	printf("between Tsys from polarizations A and B cal data.\n\n");
      } 
      else {
	printf("Will use our chosen minimum threshold value of  %5.1f%% ",
	       CalCmd->ChooseMethod);
	printf("difference between Tsys from polarizations A and B cal ");
	printf("data.\n\n");
      }
      if(CalCmd->NBadThreshP){
	printf("Will also allow only %d bad polA/polB Tsys ratios ",
	       CalCmd->NBadThresh);
	printf("(channel-to-channel, for each of ON and OFF source cals ");
	printf("separately) before reverting to constant Tsys method.\n\n");
      }      
    }
    else{
      if(CalCmd->NBadThreshP){
	printf("Cannot choose -nbad option if not using -choose option!");
	printf("Exiting...\n");     
	fflush(stdout);
	return -1;
     }
    }
  }


  /* Now get Tcal and Gain values from command line */
  for (pol=0;pol<2;pol++){
    if(CalCmd->TcalP)
      CalMode->Tcal[pol] = CalCmd->Tcal[pol];
    else
      CalMode->Tcal[pol] = -1.0;
  }
  if(CalCmd->GainP)
    CalMode->Gain = CalCmd->Gain;
  else
    CalMode->Gain = -1.0;
  if(CalCmd->GainOnPulsarP)
    CalMode->GainOnPulsar = CalCmd->GainOnPulsar;
  else
    CalMode->GainOnPulsar = CalMode->Gain;
  printf("Will assume (Gain at pulsar)/(Gain at continuum source) is: %f\n",
	CalMode->GainOnPulsar/CalMode->Gain);
  
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
       
  /* get source name */
  strcpy(RunMode->Source,CalHdr->target.PSRName); 


  /* Make sure that this *is* a cal scan! */
  if(strcmp(CalHdr->gen.ObsMode,"CAL") != 0){
    printf("File %s is NOT a cal scan according to header!  Exiting...\n",
	   RunMode->Infile);fflush(stdout);
    return -2;
  }

  /*  make sure that the Source Name has a "c" in front to distinguish */
  if (strncmp(&CalHdr->target.PSRName[0],"c",1) != 0){
    sprintf(RunMode->Source,"c%s",CalHdr->target.PSRName);    
  }
  else{
    strcpy(RunMode->Source,CalHdr->target.PSRName);fflush(stdout);
  }
  

  /* just in case */
  for(i=0;i<strlen(RunMode->Source);i++){
    if(!strncmp(&RunMode->Source[i],"_",1))
      strcpy(&RunMode->Source[i],"\0");
  }


  /* If no outfile root is given on the command line, use PSR name and some of  
   * the input file name to create one */
  LastSlashIndex = -1;
  if(!(CalCmd->OutfileRootP)) {
    for(i=0;i<strlen(RunMode->Infile);i++){
      if(!strncmp(&RunMode->Infile[i],"/",1))
	LastSlashIndex = i;
    }
    strncpy(TempChar, &RunMode->Infile[LastSlashIndex+1],12);
    strcpy(&TempChar[12],"\0");
    sprintf(RunMode->OutfileRoot,"%s.%s",RunMode->Source,TempChar);
  }
  else {
    strcpy(RunMode->OutfileRoot,CalCmd->OutfileRoot);
  }

  RunMode->FSky = CalHdr->obs.FSkyCent;

  /**** Integer Options ****/
  RunMode->NBins = CalHdr->redn.RNBinTimeDump;
  RunMode->BinDown = 0; 

  RunMode->NBinsOut = RunMode->NBins;

  /* Make sure that number of dumps is greater than zero! */
  if(RunMode->NDumps <= 0){
    printf("ERROR:  File %s has no data! Exiting\n",RunMode->Infile);fflush(stdout);
    return -3;
  } 

  /* Figure out channels to omit in getting phases */
  CalMode->NOmit=0;
  if(CalCmd->FOmitP){
    CalMode->NOmit = CalCmd->FOmitC;
    for (i=0;i<CalMode->NOmit;i++)
      CalMode->FOmit[i] = CalCmd->FOmit[i];
  }

  /* See if cal phase (where it switches from off to on (or vice versa) 
     is given by user */
  if (CalCmd->ForceP) {
    CalMode->ForcePhase=(double)CalCmd->Force;
    printf("Cal switch phase of %3.1f is chosen to be forced.\n\n",
	   CalCmd->Force);
  }
  else {
    CalMode->ForcePhase=-1.0;
    printf("Cal switch phase will be found from the data.\n\n");
  }


  return 0;

}








int GetContOpt(struct RunVars *RunMode, struct CalVars *CalMode, 
	      Cmdline *CalCmd, struct ASPHdr *CalHdr)
{



  /**** Float Options ****/

  if (!CalCmd->FluxP){
    printf("Must enter Flux of continuum source!  Exiting...\n");
    fflush(stdout);
    return -1;
  }
  CalMode->Flux = CalCmd->Flux;  /* Mandatory option if using continuum file */


  /**** Integer Options ****/
  RunMode->NBins = CalHdr->redn.RNBinTimeDump;
  RunMode->BinDown = 0; 

  RunMode->NBinsOut = RunMode->NBins;
  RunMode->OldFits = 0;


  /* Make sure that number of dumps is greater than zero!  And that it
   * agrees with what header says */
  if(RunMode->NDumps <= 0){
    printf("ERROR:  File %s has no data! Exiting\n",RunMode->Infile);fflush(stdout);
    return -3;
  } 


  /* Test command line */

  if(RunMode->Verbose) {
    printf("\nCOMMAND LINE OPTIONS:\n");
    printf("=====================\n");
    printf("Flag options:\n");
    printf("-------------\n");
    printf("RunMode->Header      = %d\n",RunMode->Header);
    printf("RunMode->FlipPA      = %d\n",RunMode->FlipPA);
    printf("RunMode->MakeRaw     = %d\n",RunMode->MakeRaw);
    printf("RunMode->Verbose     = %d\n",RunMode->Verbose);
    printf("RunMode->AddDumps    = %d\n\n",RunMode->AddDumps);
    printf("String options:\n");
    printf("---------------\n");
    printf("RunMode->Infile      = %s\n",RunMode->Infile);
    printf("RunMode->OutfileRootP = %d\n",CalCmd->OutfileRootP);
    printf("RunMode->OutfileRoot = %s\n",RunMode->OutfileRoot);
    printf("RunMode->Source      = %s\n\n",RunMode->Source);
    printf("Float options:\n");
    printf("--------------\n");
    printf("RunMode->FSky        = %f\n\n",RunMode->FSky);
    printf("Int options:\n");
    printf("------------\n");
    printf("RunMode->NBins       = %d\n",RunMode->NBins);
    printf("RunMode->NBinsOut    = %d\n\n",RunMode->NBinsOut);
    fflush(stdout);
  }

  return 0;

}
