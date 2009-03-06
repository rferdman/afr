/* Function that reads in a standard profile and rotates it to some phase */
/* Default phase is zero, but can be user-input */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


int ReadStd(struct RunVars *RunMode, struct RunVars *StdRunMode, 
	    struct StdProfs *StdProfile, double *TotWgt)
{

  int             i,jnq,bins;
  int             FoundProf;
  FILE            *fpstd;
  char            name[16],line[200], StdFileName[64], HeadLine[256];;
  double          freq,totwgtin;
  //  struct RunVars  StdRunMode;
  struct StdProfs TempProf;

  FoundProf = 0;
  RunMode->MakeStd = 0;
  /* Open standard profile file */
  if((fpstd = fopen(RunMode->Stdfile,"r")) == 0) {
    //    printf("Cannot open %s.  Will use highest signal-to-noise\n",
    //	   RunMode->Stdfile);
    // printf(" profile as standard profile...\n"); 
    RunMode->MakeStd = 1;
    memcpy(StdRunMode, RunMode, sizeof(struct RunVars));
    StdRunMode->Header = 0;
    return 1;
  }

  
  /* Read header line */
  while(fscanf(fpstd,"%lf%d%8s%lf\n",&freq,&bins,name,&totwgtin) != EOF) {
	printf("freq = %lf, bins=%d, name = %s, Source = %s\n",freq,bins,name,RunMode->Source);
 
    /* If correct Standard Prof not yet found, zero arrays again: */
    FZero(StdProfile->rstds,NBINMAX);
    FZero(StdProfile->rstdq,NBINMAX);
    FZero(StdProfile->rstdu,NBINMAX);
    FZero(StdProfile->rstdv,NBINMAX);
    FZero(StdProfile->stdlin,NBINMAX);
    FZero(StdProfile->stdphi,NBINMAX);
    FZero(StdProfile->stdphierr,NBINMAX);

    if(strncmp(RunMode->Source,"\0",1) != 0){
      if(strncmp(name,RunMode->Source,6) != 0) {
	printf("Wrong pulsar in %s.\n",RunMode->Stdfile);
	exit(1);
      }
      else {
	/* Read in profile: S, Q, U, V, L, phi, err(phi)*/
	for(i=0;i<bins;i++) {
	  fgets(line,200,fpstd);
	  sscanf(line,"%d%f%f%f%f%f%f%f\n",&jnq,&StdProfile->rstds[i],
		 &StdProfile->rstdq[i],&StdProfile->rstdu[i],
		 &StdProfile->rstdv[i],&StdProfile->stdlin[i],
		 &StdProfile->stdphi[i],&StdProfile->stdphierr[i]);
	}
      }
    }
    /* Allow for different bins in standard profile from input profile */
    /* Allow for 5% slop in frequency - enough for L-band?  */
    /* If we have the right # bins and frequency, get out of 
       this loop, we're done:*/
/*     if((bins == RunMode->NBins) && ((int)freq == (int)(RunMode->FSky))) */
    if ((fabs(freq - RunMode->FSky) < 0.05*RunMode->FSky)){
      FoundProf = 1;
      if (bins == RunMode->NBinsOut) break;
    }

  }
  fclose(fpstd);

  if(!FoundProf){
    printf("FSky = %lf,  freq = %lf\n",RunMode->FSky,freq);fflush(stdout);
    printf("Could not find standard profile at desired FSky.  Exiting...\n");
    fflush(stdout);
    return -1;
  }
  

  /* If StdProf Bins are greater than Input Prof, then BinDown Standard prof.
   * If Input Prof has more bins than standard prof, then BinDown Input Prof. */

  if(bins > RunMode->NBinsOut) {
    memcpy(&StdRunMode, RunMode, sizeof(struct RunVars));
    memcpy(&TempProf, StdProfile, sizeof(struct StdProfs));
    StdRunMode->Header = 0;
    StdRunMode->NBins = bins;
    StdRunMode->NBinsOut = RunMode->NBinsOut;
    StdRunMode->BinDown = 1;
    BinDown(StdRunMode, &TempProf, StdProfile);
    MakePol(StdRunMode, StdRunMode->NBinsOut, StdProfile);
    sprintf(StdFileName,"%s.%d.new.std",RunMode->Source,StdRunMode->NBinsOut);
    printf("Standard profile has more bins than Input Profile.\n");
    printf("Binning down standard profile and outputting file %s\n.",
	   StdFileName);fflush(stdout);
    WriteStokes(StdRunMode, StdProfile, HeadLine, StdFileName);
  }
  else if(bins < RunMode->NBinsOut) {
    RunMode->NBinsOut = bins;
    RunMode->BinDown = 1;
    printf("Input profile has more bins than standard profile.\n");
    printf("Will bin down input profile.\n");fflush(stdout);
  }




  /* allow 5% slop in frequency - enough for L-band?  */ 
  /*   if((bins != run_mode->nbins) || (fabs(freq - run_mode->fsky) > 0.05*run_mode->fsky)) { */   
  /*   printf("Std. profile of %1d bins at %5.0f MHz not found in %s\n",run_mode->nbins,run_mode->fsky,stdfile); */
  /*   exit(1); */
  /* } */

  /* Read in weight of existing standard profile */
  printf("Std. profile %s found\n",RunMode->Stdfile);
  *TotWgt = totwgtin;

  cprofc(StdProfile->rstds,RunMode->NBins,
	 StdProfile->stdamp,StdProfile->stdpha);
  cprofc(StdProfile->rstds,RunMode->NBins,
	 StdProfile->stdamps,StdProfile->stdphas);
  cprofc(StdProfile->rstdq,RunMode->NBins,
	 StdProfile->stdampq,StdProfile->stdphaq);
  cprofc(StdProfile->rstdu,RunMode->NBins,
	 StdProfile->stdampu,StdProfile->stdphau);
  cprofc(StdProfile->rstdv,RunMode->NBins,
	 StdProfile->stdampv,StdProfile->stdphav);



/*   fclose(fpstdout); */
  return 0;
}


int MakeStd(struct ASPHdr *InHdr, struct RunVars *StdRunMode, 
	    struct StdProfs *InputProfs, struct StdProfs *StdProfile, 
	    int FinalDump)
{
  static double MaxSNR=0.0;
  char          StdHead[128];

  if (InputProfs->SNR > MaxSNR){
    MaxSNR = InputProfs->SNR;
    if(StdRunMode->BinDown)
      BinDown(StdRunMode, InputProfs, StdProfile);
    else
      memcpy(StdProfile, InputProfs, sizeof(struct StdProfs));
  }
  /* if we are on the last input profile, wrap up std. prof-making */
  if(FinalDump){
    printf("Finishing std. prof-making...\n");
    printf("   Bins = %d\n",StdRunMode->NBinsOut);fflush(stdout);
    sprintf(StdHead,"%4.0f. %4d   %7s            1.000",
	    StdRunMode->FSky, StdRunMode->NBinsOut, StdRunMode->Source);
    MakePol(StdRunMode, StdRunMode->NBinsOut, StdProfile);
    WriteStokes(StdRunMode, StdProfile, StdHead, StdRunMode->Stdfile);
/*     printf("BLAH!\n");fflush(stdout); */
    cprofc(StdProfile->rstds,StdRunMode->NBinsOut,StdProfile->stdamp,
	   StdProfile->stdpha);
    cprofc(StdProfile->rstds,StdRunMode->NBinsOut,StdProfile->stdamps,
	   StdProfile->stdphas);
    cprofc(StdProfile->rstdq,StdRunMode->NBinsOut,StdProfile->stdampq,
	   StdProfile->stdphaq);
    cprofc(StdProfile->rstdu,StdRunMode->NBinsOut,StdProfile->stdampu,
	   StdProfile->stdphau);
    cprofc(StdProfile->rstdv,StdRunMode->NBinsOut,StdProfile->stdampv,
	   StdProfile->stdphav);
  }
 
  return 0;
  
}
