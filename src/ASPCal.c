/********************************************************************/
/*                                                                  */
/*  ASPCal:                                                         */
/*                                                                  */
/*  Will perform calibrations on cal scan data in order to obtain:  */
/*                                                                  */
/*     - Gain calibration factors through Tcal or continuum source  */
/*       observations                                               */
/*     - Calculations of Tsys, Tant                                 */
/*                                                                  */
/* Robert Ferdman, University of British Columbia, August 24, 2004  */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "fitsio.h"
#include "ASPCommon.h"
#include "CalCmdLine.h"
#include "ASPCal.h"


int main(int argc, char **argv)
{
  
  int            i, j, i_loop, n_loop, pol, status=0, NumHDU, LastSlashIndex;

  double         *ASquared[2*NCHMAX], *BSquared[2*NCHMAX];
  double         *ReAconjB[2*NCHMAX], *ImAconjB[2*NCHMAX];

  struct ASPHdr  CalHdr, ContHdr[2];
  struct SubHdr  *CalSubHdr, *ContSubHdr[2];
  struct RunVars CalRunMode, ContRunMode[2];
  struct CalVars CalMode;


  //  int            PhaseBin[2], ContPhaseBin[2][2];
  int            OnBin[2],OffBin[2], ContOnBin[2][2],ContOffBin[2][2];

  int            ONSource=0, OFFSource=0; // just to initialize for compilation
  int            NumChansFound;
  int            NBadRatio=0;
  int            SkipChan[NCHMAX];
  //  double         HighAvg[2], LowAvg[2];
  double         OnAvg[2], OffAvg[2];
  double         CalHeight[NCHMAX][2], ContCalHeight[2], ContCalFrac[2];
  // double         DeltaOff, DeltaOn;

  double         JyPerCount[NCHMAX][4], JyPerCal[NCHMAX][4];
  double         Tsys[2][2], Tant, Tcal;
  double         TsysRatio[NCHMAX][2], AvgTsysRatio[2];
  // double         Gain;
  // double         ThetaBB;

  char           ProgName[256];
  char           OutfileRoot[256];
  char           TempChar[20];
  Cmdline        *CalCmd;

  fitsfile       *Fcal, *Fcont[2];
  char           CalOutFile[64];
  FILE           *Fcalout;

  int            Diagnose=0;
  double         JyPerCountCont[NCHMAX][4];

  FILE           *Fcheck, *Ftsys;

  /* Get command line variables */
  CalCmd = parseCmdline(argc, argv);
  /*   showOptionValues();  */

  /* Normally use this somewhere, and not showOptionValues */
  CalCmd->tool = CalCmd->tool;
                                                                             
  /* Store program name */
  strcpy(ProgName,argv[0]);  

  /* Zero out SkipChan array */
  IZero(&SkipChan[0], NCHMAX);

  /*** PUT IN READING OF PULSAR CAL FILE HERE ***/

  if (CalCmd->CalfileP){

    /* Initialize header variables */
    InitPars(&CalHdr);

    /* Grab pulsar cal file name */
    strcpy(CalRunMode.Infile,CalCmd->Calfile);

    /* Open pulsar cal fits file */
    status=0;

    if(fits_open_file(&Fcal, CalRunMode.Infile, READONLY, &status)){
      printf("Error opening FITS cal file %s !!!\n", CalRunMode.Infile);
      exit(1);
    }

    /* Get number of HDUs in fits file */
    fits_get_num_hdus(Fcal, &NumHDU, &status);
    //    CalRunMode.NDumps = NumHDU-3;  

    /* Read in values for header variables */
    if(ReadASPHdr(&CalHdr, Fcal) < 0){
      printf("Unable to read Header from CAL file %s.  Exiting...\n",
	     CalRunMode.Infile);
      exit(2);
    };

    if(!strcmp(CalHdr.gen.HdrVer,"Ver1.0")){
      CalRunMode.NDumps += NumHDU-3;  /* the "3" is temporary, depending on how 
				     many non-data tables we will be using */
    }
    else if(!strcmp(CalHdr.gen.HdrVer,"Ver1.0.1")){
      CalRunMode.NDumps += (NumHDU-3)/2;
    }
    else{
      printf("Do not recognize FITS file version number in header.\n");
      printf("This header %s. Exiting...\n",CalHdr.gen.HdrVer);fflush(stdout);
      exit(3);
    }

    printf("\n==========================\n");
    printf("ASP FITS Header %s\n",CalHdr.gen.HdrVer);
    printf("==========================\n\n");fflush(stdout);
    
    printf("Input cal file for PSR %s:\n    %s\n\n",
	   CalHdr.target.PSRName,CalRunMode.Infile);fflush(stdout);
    
    printf("Centre Frequency: %6.1lf MHz\n\n",CalHdr.obs.FSkyCent);fflush(stdout);
    


    /* Get options */
    if(GetCalOpt(&CalRunMode, &CalMode, CalCmd, &CalHdr) < 0){
      printf("Unable to get options.  Exiting...\n");
      exit(3);
    }

    /* Read in cal data */
    if(GetCalData(&CalHdr, CalSubHdr, &CalRunMode, Fcal, 
		  &ASquared[0], &BSquared[0], &ReAconjB[0], &ImAconjB[0]) < 0) {
      printf("Error: Could not read cal data from file %s.  Exiting...\n",
	     CalRunMode.Infile);fflush(stdout);
      exit(4);
    }




    if(CalCmd->OutfileRootP){
      sprintf(CalOutFile,"%s.cal",CalCmd->OutfileRoot);
    }
    else{
      LastSlashIndex = -1;
      for(i=0;i<strlen(CalRunMode.Infile);i++){
	if(!strncmp(&CalRunMode.Infile[i],"/",1))
	  LastSlashIndex = i;
      }
      strncpy(TempChar, &CalRunMode.Infile[LastSlashIndex+1],12);
      strcpy(&TempChar[12],"\0");
     
      sprintf(OutfileRoot,"%s.%s.%s",&CalRunMode.Source[1],TempChar,CalHdr.obs.ObsvtyCode);
      if(CalRunMode.Verbose)
	printf("Infile is %s, OutFileRoot is %s\n",CalRunMode.Infile, 
	       OutfileRoot);fflush(stdout);      
      
      strcpy(CalOutFile,"\0");
      sprintf(CalOutFile,"%s.cal",OutfileRoot);
      
      //    sprintf(CalOutFile,"%s.%5.5d.%.d.cal",CalRunMode.Source, CalHdr.obs.IMJDStart, (int)CalHdr.obs.FSkyCent);     
    }
    if((Fcalout = fopen(CalOutFile,"w")) == NULL){
      printf("Could not open file %s.  Exiting...\n",CalOutFile);fflush(stdout);
      exit(6);
    }

 
    /* Find ranges of bins for on and off for each channel -- default is all */
    if (GetPhases(&CalHdr, &CalRunMode, &CalMode, 
		  &ASquared[0], &BSquared[0], &OnBin[0], &OffBin[0]) <0){
      printf("Problem getting on/off phases for cal file %s. Exiting...\n",
	     CalRunMode.Infile);fflush(stdout);
      exit(5);
    }


    if(CalRunMode.Verbose)
      printf("PULSAR CAL SCAN FILE:\n\n");
    /* Open output file which will be read into ASPFitsReader */

    for(i=0;i<CalHdr.obs.NChan;i++){

    
    /* Open output file which will be read into ASPFitsReader */
      if(CalRunMode.Verbose)
	printf("CHAN %d = %lf MHz:\n",i,CalHdr.obs.ChanFreq[i]);

      for(pol=0;pol<2;pol++){
      

	/* Find average of cal on and off separately */
	if(pol==0)
	  CalHeight[i][pol] = GetCalHeight(ASquared[i],CalRunMode.NBins,OnBin,OffBin,
					   &OnAvg[0],&OffAvg[0]);
	else
	  CalHeight[i][pol] = GetCalHeight(BSquared[i],CalRunMode.NBins,OnBin,OffBin,
					   &OnAvg[0],&OffAvg[0]);
      
	/* If cal is weird and gives negative values, then skip in and don't 
	   write JyPerCount to file.  AFR will omit that channel 
	   automatically.*/
	if (CalHeight[i][pol] < 0.) {
	  printf("WARNING: Did not calculate CalHeight correctly. ");
	  printf("Skipping and will not write channel %.1lf to file.\n",
		 CalHdr.obs.ChanFreq[i]);
	  SkipChan[i]=1;
	  // return -9;
	}
	else {
	  if(CalRunMode.Verbose)
	    printf("OnAvg[%d] = %lf, OffAvg[%d] = %lf\n",
		   pol,OnAvg[0],pol,OffAvg[0]);
	  
	  /* Calculate only if we are not using continuum cals */
	  if (!CalCmd->ContfileP)
	    JyPerCount[i][pol] = CalMode.Tcal[pol]/(CalMode.Gain*CalHeight[i][pol]); 
	}
      }
      /* Calculate Jy/count for polarization cross-products, again, 
       * if not using continuum cal file */
      if (!CalCmd->ContfileP){
	/* Calculate and print to file only if channel is good and not being
	   skipped */
	if (!SkipChan[i]){
	JyPerCount[i][2] = JyPerCount[i][3] 
	  = sqrt(JyPerCount[i][0]*JyPerCount[i][1]);          

	if(CalRunMode.Verbose)
	  printf("JyPerCal: %lf  %lf  %lf  %lf\n\n",
		 JyPerCount[i][0],JyPerCount[i][1],
		 JyPerCount[i][2],JyPerCount[i][3]);
      
	/* write cal factors to file */
	fprintf(Fcalout,"%lf  %lf  %lf  %lf  %lf\n",CalHdr.obs.ChanFreq[i],
		JyPerCount[i][0],JyPerCount[i][1],
		JyPerCount[i][2],JyPerCount[i][3] );
	}
      }
    }

    fits_close_file(Fcal, &status);

  }
  


  /*** END PULSAR CAL FILE SECTION ***/
  
  
  
  /*** BEGIN READING AND PROCESSING CONTINUUM CAL FILE ***/
  
  /* Read in each cal file (ON and OFF source) separately */
  
  if (CalCmd->ContfileP) {
 
    for (i=0;i<2;i++){ // ON and OFF files

      /* Initialize header variables */
      InitPars(&ContHdr[i]);

      /* Grab continuum cal file name */
      strcpy(ContRunMode[i].Infile,"\0");
      strcpy(ContRunMode[i].Infile,CalCmd->Contfile[i]);
     
      /* Open continuum cal fits file */
      status=0;

      if(fits_open_file(&Fcont[i], ContRunMode[i].Infile, READONLY, &status)){
	printf("Error opening FITS cal file %s !!!\n", ContRunMode[i].Infile);
	exit(1);
      }
      /* Get number of HDUs in fits file */
      fits_get_num_hdus(Fcont[i], &NumHDU, &status);
     /*ContRunMode[i].NDumps = NumHDU-3;*//* the "3" is temporary, depending on 
					  * how many non-data tables we will be
					  * using */

      /* Read in values for header variables */
      if(ReadASPHdr(&ContHdr[i], Fcont[i]) < 0){
	printf("Unable to read Header from Cont. CAL file %s.  Exiting...\n",
	       ContRunMode[i].Infile);
	exit(2);
      };
    if(!strcmp(CalHdr.gen.HdrVer,"Ver1.0")){
      ContRunMode[i].NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				     many non-data tables we will be using */
    }
    else if(!strcmp(CalHdr.gen.HdrVer,"Ver1.0.1")){
      ContRunMode[i].NDumps = (NumHDU-3)/2;
    }
    else{
      printf("Do not recognize FITS file version number in header.\n");
      printf("This header %s. Exiting...\n",CalHdr.gen.HdrVer);fflush(stdout);
      exit(3);
    }


    /* If second of two files, make sure that the continuum cal files are for 
       the same source */
    
    /* Also make sure that all channels in continuum files are represented 
       in pulsar cal file (== vice-versa) */
    

    /* Get options */
    if(GetContOpt(&ContRunMode[i], &CalMode, CalCmd, &ContHdr[i]) < 0){
      printf("Unable to get options.  Exiting...\n");
      exit(3);
    }
    
    /* Read in cal files (ON and OFF cal source) */
    
    if(GetCalData(&ContHdr[i], ContSubHdr[i], &ContRunMode[i], Fcont[i], 
		  &ASquared[i*NCHMAX], &BSquared[i*NCHMAX],
		  &ReAconjB[i*NCHMAX], &ImAconjB[i*NCHMAX]) < 0) {
      printf("Error: Could not read cal data from file %s.  Exiting...\n",
	     ContRunMode[i].Infile);fflush(stdout);
      exit(4);
    }
    
    /* Add all data from all channels together in order to find phases 
     * at which cal in on and off, and skip omitted channels */
    if (GetPhases(&ContHdr[i], &ContRunMode[i], &CalMode, 
		  &ASquared[i*NCHMAX], &BSquared[i*NCHMAX], 
		  ContOnBin[i], ContOffBin[i]) < 0){
      printf("Problem obtaining on/off phases for cal file %s. Exiting...\n",
	     ContRunMode[i].Infile);fflush(stdout);
      exit(5);
    }
    
    if(fits_close_file(Fcont[i], &status)){
      printf("Problems closing fits file %s.  Exiting...\n",
	     ContRunMode[i].Infile);fflush(stdout);
      exit(7);
    }
    
    }
    
    /*** END READING AND PROCESSING CONTINUUM CAL FILE ***/
    
    /* Now I have cal profiles for each of ON and OFF continuum source */
    
    /* Open output cal file */
    
    
    /*** BEGIN CONTINUUM CALIBRATION ***/
    
    if(CalRunMode.Verbose){
      Fcheck = fopen("check_cont.dat", "w");
      Ftsys =  fopen("tsys_ratio.dat", "w");
    }

    if(CalCmd->ConstTsysP){
      printf("We are assuming CONSTANT TSYS.  Values of Tsys, Tant, Tcal ");
      printf("are given in K, based on input command line Gain value.\n\n");
      printf("JyPerCount is based on the off-source cal height.  ");
      printf("It is recalculated for pulsar cal based on JyPerCal ");
      printf("found.\n\n");
    }
    else{
      if(!CalCmd->ChooseMethodP){
	printf("We are assuming CONSTANT GAIN.  Values of Tsys are based ");
	printf("on the difference in on-vs-off-source cal-off values and ");
	printf("input gain (as are Tcal and Tant).\n\n");
	printf("JyPerCount is based on the off-source cal height.  ");
	printf("It is recalculated for pulsar cal based on JyPerCal ");
	printf("found.\n\n");
      }
    }

    /* Find Scal, ie the Jy/cal using the continuum source */

    printf("\n\nResults from continuum source calibrator:\n");
    printf("-----------------------------------------\n");

    /* Prepare to go through loop twice if letting ASPCal choose cal method */ 
    if(CalCmd->ChooseMethodP){
      /* Zero out TsysRatio averages if choosing cal method */
      DZero(AvgTsysRatio,2);
      n_loop=2;
    }
    else{
      n_loop=1;
    }

    for(i_loop=0;i_loop<n_loop;i_loop++){

    for (pol=0;pol<2;pol++) {
    
      if (pol == 0){ 
	printf("\nPolarization A:\n\n");
	if(CalRunMode.Verbose)
	  fprintf(Fcheck,"\nPolarization A:\n\n");
      }
      else {
	printf("\nPolarization B:\n\n");
	if(CalRunMode.Verbose)
	  fprintf(Fcheck,"\nPolarization B:\n\n");
      }

      printf(" Freq(MHz)  JyPerCal   JyPerCount   TsysON(K)  TsysOFF(K)  Tant(K)  Tcal(K)\n");
      printf(" ---------  --------   ----------   ---------  ----------  -------  -------\n");


      for (j=0;j<ContHdr[0].obs.NChan;j++){

	if(CalRunMode.Verbose)
	  fprintf(Fcheck,"CHAN %d = %lf MHz:\n",j,ContHdr[0].obs.ChanFreq[j]);


	for (i=0;i<2;i++) {

	  if (pol == 0) 
	    ContCalHeight[i] = GetCalHeight(ASquared[i*NCHMAX + j], 
					    ContRunMode[i].NBins,
					    ContOnBin[i],ContOffBin[i],
					    &OnAvg[i], &OffAvg[i]);
	  else
	    ContCalHeight[i] = GetCalHeight(BSquared[i*NCHMAX + j], 
					    ContRunMode[i].NBins,
					    ContOnBin[i],ContOffBin[i],
					    &OnAvg[i], &OffAvg[i]);
	  if (ContCalHeight[i] < 0.) {
	    printf("WARNING: Did not calculate CalHeight correctly for ");
	    printf("continuum cal scan %s. ",ContRunMode[i].Infile);
	    printf("Skipping and will not write channel %.1lf to file.\n",
		 ContHdr[0].obs.ChanFreq[j]);
	    SkipChan[j]=1;
	  }

	/* Calculte calheight/baseline ratio for each scan and use that to 
	   determine ON and OFF Source */
	  ContCalFrac[i] = ContCalHeight[i]/OffAvg[i];

	}

	

	//	if(OffAvg[0] > OffAvg[1]){
	if(ContCalFrac[0] < ContCalFrac[1]){
	  if((j>0 && pol==0) && ONSource==1) { // i.e. opposite to last channel result -- also both pol's must agree
	    printf("Ambiguity as to which continuum scan is ON source and ");
	    printf("which is OFF source. Check command line and try again. ");
	    printf("Exiting...\n");
	    exit(2);
	  }
	  ONSource = 0; 
	  OFFSource = 1;
	}
	else{
	  if((j>0 && pol==0) && ONSource==0) { // i.e. opposite to last channel result -- also both pol's must agree
	    printf("Ambiguity as to which continuum scan is ON source and ");
	    printf("which is OFF source. Check command line and try again. ");
	    printf("Exiting...\n");
	    exit(2);
	  }
	  ONSource = 1; 
	  OFFSource = 0;
	}

	if(CalRunMode.Verbose){
	  fprintf(Fcheck,"OnAvg[ONSource] =  %lf, OffAvg[ONSource] =  %lf\n",
		  OnAvg[ONSource], OffAvg[ONSource]);
	  fprintf(Fcheck,"OnAvg[OFFSource] = %lf, OffAvg[OFFSource] = %lf\n",
		  OnAvg[OFFSource], OffAvg[OFFSource]);
	}

	if(CalCmd->ConstTsysP) {
	  if (GetCalTsys(&CalMode,
			 OnAvg, OffAvg, ContCalHeight, 
			 ONSource, OFFSource, 
			 &JyPerCal[j][pol], &JyPerCount[j][pol], 
			 Tsys[pol], &Tant, &Tcal) < 0){
	    printf("Could not get Cal Temperatures through the constant gain ");
	    printf("assumption.  Exiting...");
	    fflush(stdout);
	    exit(2);
	  }
	}
	else{
	  if (GetCalGain(&CalMode,
			 OnAvg, OffAvg, ContCalHeight, 
			 ONSource, OFFSource, 
			 &JyPerCal[j][pol], &JyPerCount[j][pol], 
			 Tsys[pol], &Tant, &Tcal) < 0){
	    printf("Could not get Cal Temperatures through the constant gain ");
	    printf("assumption.  Exiting...");
	    fflush(stdout);
	    exit(2);
	  }

	  /* If we want ASPCal to choose the method then keep a record of Tsys
	     ratios */
	}	

	printf("%10.1lf  %8.2lf  %11.2lf  %10.3lf  %10.3lf  %7.1lf  %7.1lf\n",
	       ContHdr[0].obs.ChanFreq[j],JyPerCal[j][pol], JyPerCount[j][pol],
	       Tsys[pol][ONSource],Tsys[pol][OFFSource],Tant,Tcal);
	fflush(stdout);

	/* Calculate Tsys Ratio if we want ASPCal to choose cal method */
	if(CalCmd->ChooseMethodP){
	  if (pol==0){ // use TsysRatio array to store Tsys's for pol==0:
	    TsysRatio[j][ONSource] = Tsys[0][ONSource];
	    TsysRatio[j][OFFSource] = Tsys[0][OFFSource];
	  }
	  else{ // i.e. pol==1 -- overwrite, now will be actual TsysRatio:
	    TsysRatio[j][ONSource] = 100.* fabs( (Tsys[1][ONSource]/TsysRatio[j][ONSource]) - 1. );
	    TsysRatio[j][OFFSource] = 100.* fabs( (Tsys[1][OFFSource]/TsysRatio[j][OFFSource]) - 1. );
	    if(TsysRatio[j][ONSource] > CalCmd->ChooseMethod ||
	       TsysRatio[j][OFFSource] > CalCmd->ChooseMethod){
	      NBadRatio++;
	    }
	    if(CalRunMode.Verbose){
	      fprintf(Ftsys,"%6.1lf   %lf   %lf  %d  %d\n", ContHdr[0].obs.ChanFreq[j], TsysRatio[j][ONSource], TsysRatio[j][OFFSource], ONSource, OFFSource);
	    }

	    /* Keep running total of ratios */
	    AvgTsysRatio[ONSource] += TsysRatio[j][ONSource];
	    AvgTsysRatio[OFFSource] += TsysRatio[j][OFFSource];
	  }

	}
	
      }
    }

    /* Now test mean Tsys ratio for each polarization.  If it exceeds 
       threshold, then go through loop again and do constant Tsys method.  If 
       is within threshold, keep constant-gain results, break out of loop,  
       and write out cal file to feed into ASPFitsReader */
    if(CalCmd->ChooseMethodP){
      AvgTsysRatio[ONSource] /= (double)ContHdr[0].obs.NChan;
      AvgTsysRatio[OFFSource] /= (double)ContHdr[0].obs.NChan;
      if(AvgTsysRatio[ONSource] > CalCmd->ChooseMethod ||
	 AvgTsysRatio[OFFSource] > CalCmd->ChooseMethod){
	printf("\nMean Tsys ratio for ON Source cal:  %5.1lf%%\n",
	       AvgTsysRatio[ONSource]);
	printf("Mean Tsys ratio for OFF Source cal: %5.1lf%%\n\n",
	       AvgTsysRatio[OFFSource]);
	printf("At least one of these is larger than the %5.1lf%% threshold.  ",
	       CalCmd->ChooseMethod);
	printf("Reverting to the constant Tsys method...\n\n");
	sleep(4);
	CalCmd->ConstTsysP=1;
	CalCmd->ChooseMethodP=0;
      }
      else if (CalCmd->NBadThreshP && (NBadRatio > CalCmd->NBadThresh) ) {
	printf("\nThere are at least %d instances of bad Tsys ratio(s), ",
	       NBadRatio);
	printf("which exceed(s) the %5.1lf%% threshold.  This is larger ",
	       CalCmd->ChooseMethod);
	printf("than the allowable amount of %d bad ratios specified on the ",
	       CalCmd->NBadThresh);
	printf("command line. Reverting to the constant Tsys method...\n\n");
	sleep(4);
 	CalCmd->ConstTsysP=1;
	CalCmd->ChooseMethodP=0;
      }
      else{
	printf("\nMean Tsys ratio for ON Source cal:  %5.1lf%%\n",
	       AvgTsysRatio[ONSource]);
	printf("Mean Tsys ratio for OFF Source cal: %5.1lf%%\n\n",
	       AvgTsysRatio[OFFSource]);
	printf("Both are within the %5.1lf%% threshold. Will continue ",
	       CalCmd->ChooseMethod);
	printf("using constant gain assumption.\n\n");
	break;
      }
    }
    }

    if(CalRunMode.Verbose){
      fclose(Fcheck);
      if(CalCmd->ChooseMethodP) 
	fclose(Ftsys);
    }

    printf("\nON Source %s continuum cal file:  %s\n",
	   ContHdr[ONSource].target.PSRName,ContRunMode[ONSource].Infile);
    printf("OFF Source %s continuum cal file:  %s\n\n",
	   ContHdr[OFFSource].target.PSRName,ContRunMode[OFFSource].Infile);

    
    /* Now actually calculate JyPerCount if user has provided pulsar cal file.
     * Otherwise, give warning that value calculated is only based on 
     * continuum cal scan */

    /* First make sure that we match channels, just in case we want 
     * (for whatever reason) to use a continuum cal scan with different 
     * numbers of channels */

    printf("Writing file %s\n\n",CalOutFile);fflush(stdout);

    if(CalCmd->CalfileP){ // here, output cal file is already opened...
      NumChansFound = 0;
      for(i=0;i<ContHdr[0].obs.NChan;i++){
	for (j=0;j<CalHdr.obs.NChan;j++){
	  if(CalHdr.obs.ChanFreq[j] == ContHdr[0].obs.ChanFreq[i]) {
	    NumChansFound++;
	    
	    /* Write to file if this channel is not bad data */
	    if (!SkipChan[j]) {
	      for(pol =0;pol<2;pol++){
		JyPerCount[j][pol] = JyPerCal[i][pol]/CalHeight[j][pol];
	      }
	      
	      /* Write JyPerCount's to file as we find matching frequencies */
	      JyPerCount[j][2] = JyPerCount[j][3] 
		= sqrt(JyPerCount[j][0]*JyPerCount[j][1]);          
	      fprintf(Fcalout,"%lf  %lf  %lf  %lf  %lf\n",
		      CalHdr.obs.ChanFreq[j],
		      JyPerCount[j][0],JyPerCount[j][1],
		      JyPerCount[j][2],JyPerCount[j][3] );
	    }
	    
	    break; // move on to next continuum cal frequency to match up
	  }
	}
      }        

    
      if(NumChansFound != CalHdr.obs.NChan){
	printf("Warning:  All channels in input pulsar cal data file are\n"); 
	printf("  NOT accounted for in continuum cal file. Writing only\n");
	printf("  matching frequencies to file\n\n");
      }
    }
    else{

      printf("Warning: Jy-Per-Count values are taken directly from\n");
      printf("  Continuum cal files");fflush(stdout);

      /* Open output file which will be read into ASPFitsReader */
      sprintf(CalOutFile,"%s.%5.5d.cal",ContRunMode[0].Source, 
	      ContHdr[0].obs.IMJDStart);
      if((Fcalout = fopen(CalOutFile,"w")) == NULL){
	printf("Could not open file %s.  Exiting...\n",CalOutFile);
	fflush(stdout);
	exit(6);
      }      

      /* write to output file */
      for (j=0;j<ContHdr[0].obs.NChan;j++){
	/* Write to file if this channel is not bad data */
	if (!SkipChan[j]) {
	  JyPerCount[j][2] = JyPerCount[j][3] 
	    = sqrt(JyPerCount[j][0]*JyPerCount[j][1]);          
	  fprintf(Fcalout,"%lf  %lf  %lf  %lf  %lf\n",
		  ContHdr[0].obs.ChanFreq[j],
		  JyPerCount[j][0],JyPerCount[j][1],
		  JyPerCount[j][2],JyPerCount[j][3] );
	}
      }
    }

    
   

    /*** END CONTINUUM CALIBRATION ***/

  }

  /* If error-checking, print separate file for each of ON source, OFF source, 
     and cal on and cal off, that contains A, B, A*CalFac(A), B*CalFac(B) */
  if (CalCmd->ContfileP && Diagnose){

    int chan, nfile;
    int AvgOnBin, AvgOffBin;
    double Atest_on, Atest_off, Btest_on, Btest_off, ACal, BCal;
    FILE *fptest[4];

    if((fptest[0] = fopen("test_rawAB_OFFSource_CalON.dat","w")) == NULL){
      printf("Could not open diagnostic file test_rawAB_OFFSource_CalON.dat.  Exiting...");fflush(stdout);
      exit(2);
    }
    if((fptest[1] = fopen("test_rawAB_OFFSource_CalOFF.dat","w")) == NULL){
      printf("Could not open diagnostic file test_rawAB_OFFSource_CalOFF.dat.  Exiting...");fflush(stdout);
      exit(2);
    }
    if((fptest[2] = fopen("test_rawAB_ONSource_CalON.dat","w")) == NULL){
      printf("Could not open diagnostic file test_rawAB_ONSource_CalON.dat.  Exiting...");fflush(stdout);
      exit(2);
    }
    if((fptest[3] = fopen("test_rawAB_ONSource_CalOFF.dat","w")) == NULL){
      printf("Could not open diagnostic file test_rawAB_ONSource_CalOFF.dat.  Exiting...");fflush(stdout);
      exit(2);
    }


    if(ContOnBin[OFFSource][0]<ContOnBin[OFFSource][1]){
      AvgOnBin = (int)(((float)(ContOnBin[OFFSource][0])+
			(float)(ContOnBin[OFFSource][1]))/2.0);
      if((AvgOffBin = AvgOnBin+(int)((float)ContRunMode[OFFSource].NBins/2.0))
	 > ContRunMode[OFFSource].NBins) 
	AvgOffBin = AvgOnBin-(int)((float)ContRunMode[OFFSource].NBins/2.0);

    }
    else{
      AvgOffBin = (int)(((float)(ContOffBin[OFFSource][0])+
			 (float)(ContOffBin[OFFSource][1]))/2.0);
      if((AvgOnBin = AvgOffBin+(int)((float)ContRunMode[OFFSource].NBins/2.0))
	 > ContRunMode[OFFSource].NBins) 
	AvgOnBin = AvgOffBin-(int)((float)ContRunMode[OFFSource].NBins/2.0);
    }

    printf("\nAvgONBin = %d, AvgOffBin = %d\n\n",
	   AvgOnBin,AvgOffBin);fflush(stdout);
  
    for (chan=0;chan<ContHdr[0].obs.NChan;chan++){
        /* OFF SOURCE */

	Atest_on = ASquared[OFFSource*NCHMAX + chan][AvgOnBin];
	Btest_on = BSquared[OFFSource*NCHMAX + chan][AvgOnBin];
	Atest_off = ASquared[OFFSource*NCHMAX + chan][AvgOffBin];
	Btest_off = BSquared[OFFSource*NCHMAX + chan][AvgOffBin];

	ACal = Atest_on*JyPerCountCont[chan][0];
	BCal = Btest_on*JyPerCountCont[chan][1];
	fprintf(fptest[0],"%6.1lf MHz:   %10.4lf  %10.4lf  %10.4lf  %10.4lf\n",
		ContHdr[0].obs.ChanFreq[chan], 
		Atest_on, Btest_on, ACal, BCal );

	ACal = Atest_off*JyPerCountCont[chan][0];
	BCal = Btest_off*JyPerCountCont[chan][1];
	fprintf(fptest[1],"%6.1lf MHz:   %10.4lf  %10.4lf  %10.4lf  %10.4lf\n",
		ContHdr[0].obs.ChanFreq[chan], 
		Atest_off, Btest_off, ACal, BCal );
	
	
	/* ON SOURCE */

	Atest_on = ASquared[ONSource*NCHMAX + chan][AvgOnBin];
	Btest_on = BSquared[ONSource*NCHMAX + chan][AvgOnBin];
	Atest_off = ASquared[ONSource*NCHMAX + chan][AvgOffBin];
	Btest_off = BSquared[ONSource*NCHMAX + chan][AvgOffBin];

	ACal = Atest_on*JyPerCountCont[chan][0];
	BCal = Btest_on*JyPerCountCont[chan][1];
	fprintf(fptest[2],"%6.1lf MHz:   %10.4lf  %10.4lf  %10.4lf  %10.4lf\n",
		ContHdr[0].obs.ChanFreq[chan], Atest_on, Btest_on, ACal, BCal );

	ACal = Atest_off*JyPerCountCont[chan][0];
	BCal = Btest_off*JyPerCountCont[chan][1];
	fprintf(fptest[3],"%6.1lf MHz:   %10.4lf  %10.4lf  %10.4lf  %10.4lf\n",
		ContHdr[0].obs.ChanFreq[chan], Atest_off, Btest_off, ACal, BCal );

    }
    
    for (nfile=0;nfile<4;nfile++) fclose(fptest[nfile]);

  }
  
  /* Calculate JyPerCount and make output file for ASPFitsReader to use, 
   * containing JyPerCount */

  
   fclose(Fcalout);

   printf("Completed successfully.\n\n");fflush(stdout);

  exit(0);

}


/* Find Cal Height and counts in Jy assuming constant gain for on/off source */
int GetCalGain (struct CalVars *CalMode, 
		double *OnAvg, double *OffAvg, double *ContCalHeight, 
		int ONSource, int OFFSource, 
		double *JyPerCal, double *JyPerCount, 
		double *Tsys, double *Tant, double *Tcal)
{

  double DeltaOff, DeltaOn, KPerCount;
  
  DeltaOff =  OffAvg[ONSource]  - OffAvg[OFFSource];
  DeltaOn = OnAvg[ONSource] - OnAvg[OFFSource];

  *JyPerCount = CalMode->Flux/DeltaOff;
  //	if(Diagnose) JyPerCountCont[j][pol] = JyPerCount[j][pol];
  *JyPerCal = *JyPerCount * ContCalHeight[OFFSource];

  /* Now use that to calculate Tsys -- if Tcal or Gain is given, 
   * can compare and verify */

  /* Diagnostics if Tcal and/or Gain are given on command line */
  /*  if(CalMode->Tcal[pol] > 0.0){
    KPerCount  = (CalMode.Tcal[pol]/ContCalHeight[OFFSource]);
  
    Tsys[ONSource]  = KPerCount * OffAvg[ONSource];
    Tsys[OFFSource] = KPerCount * OffAvg[OFFSource];

    Tant = KPerCount * DeltaOff;
    Gain = Tant/CalMode.Flux;
    }  */

  /* Telescope gain is a mandatory command line argument so we can
     calculate this */
  *Tant = CalMode->Gain * CalMode->Flux;
  KPerCount = *Tant/DeltaOff;
	
  Tsys[ONSource]  = KPerCount * OffAvg[ONSource];
  Tsys[OFFSource] = KPerCount * OffAvg[OFFSource];
	
  *Tcal = KPerCount * ContCalHeight[OFFSource];
	
  return 0;

}


/* Find Cal Height and counts in Jy assuming constant Tsys for on/off source */
/* Based on Willem's psrchive routine, fluxcal */
int GetCalTsys(struct CalVars *CalMode, 
		double *OnAvg, double *OffAvg, double *ContCalHeight, 
		int ONSource, int OFFSource, 
		double *JyPerCal, double *JyPerCount, 
		double *Tsys, double *Tant, double *Tcal)
{

  double fon, foff, ffactor;
  
  fon = (OnAvg[ONSource]/OffAvg[ONSource]) - 1.;
  foff = (OnAvg[OFFSource]/OffAvg[OFFSource]) - 1.;

  ffactor = (1./fon) - (1./foff);

  /*  printf("\n\n FON = %lf,   FOFF = %lf,   FFACTOR = %lf\n\n", 
      fon, foff, ffactor);
      fflush(stdout);  */
  
  *JyPerCal = CalMode->Flux / ffactor;
  
  /* Solve for Tsys and convert to K based on estimated gain given on command 
     line */
  Tsys[ONSource] = Tsys[OFFSource] = CalMode->Gain * (*JyPerCal / foff);

  /* Jy per Count based on height of off-source cal height  -- should be 
     constant according to this method? JPerCount is calculated separately 
     for pulsar cal anyhow */
  *JyPerCount  = (*JyPerCal)/ContCalHeight[OFFSource];

  /* Tant based on a constant estimated gain from command line */
  *Tant = CalMode->Gain * CalMode->Flux;

  /* Tcal based on estimated constant gain from command line */
  *Tcal = *JyPerCal * CalMode->Gain;

  return 0;
}
