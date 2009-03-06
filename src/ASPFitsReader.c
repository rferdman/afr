/* ========================================================================= */
/* ASPFitsReader                                                             */
/*    Reads in fits output from ASP, performs calibration, and outputs       */
/*    a calibrated Stokes parameter FITS file                                */
/*                                                                           */
/* R. Ferdman, 18-January-2004, University of British Columbia.              */
/*             Latest Version -- Dec. 13, 2004                               */
/* ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "ASPCommon.h"
#include "ASPFitsReader.h"
#include "polyco.h"

int main(int argc, char **argv)
{
  int             i,i_dump_in,i_dump_out,i_chan_in,i_chan_out,i_out,i_bin; 
  int             MinAddDump, MaxAddDump, OutChanIndex, denom;
  int             SubOutDumpIndex, SubOutChanIndex;
  int             hdutype, status=0, q;
  int             RootIndex=0, NewDump=0, NumDumps=-1;
  int             NFirstTable, NumHDU, FileNo, NFilesOut, FileOutNo;
  long            NPtsProf=0;
  char            ProgName[256], *HeadLine[64];
  char            *AddChansHead[64], *OutputHead[64];
  char            FitsFile[128], FitsFileOut[128];
  fitsfile        *Fin, **Fout;
  int             fitsstatus=0;
  struct ASPHdr   InHdr, OutHdr;
  struct SubHdr   *SubInHdr, SubOutHdr; //, *SubDumpHdr, 
  struct RunVars  RunMode, StdRunMode;
  struct CalVars  CalMode;
  struct StdProfs StdProfile, *InputProfs, TempOutProfs, StokesProfs;
  struct StdProfs *AddChansProfs, *OutputProfs; 
  double          **ASquared, **BSquared, **ReAconjB, **ImAconjB;
  int             **SampleCount;
  double          *JyPerCount[NCHMAX];
  double          TotWgt;
  char            Stokesfile[256], TempInFile[256];

  double          StartMJD;
  struct Polyco   *Polycos;
  int             n_poly;

  Cmdline         *Cmd;


  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);
  /*  showOptionValues();  */

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;
                                                                             
  /* Store program name */
  strcpy(ProgName,argv[0]);  

  /* Initialize header variables */
  InitPars(&InHdr);

  /* Grab infile name */
  strcpy(RunMode.Infile,Cmd->Infile); 

  /* Malloc fitsfile descriptor to number of data files */


  /* Open fits data file */

  /* Check that file ends in ".asp" */
  for (i=0;i<strlen(RunMode.Infile);i++){
    if(!(strcmp(&RunMode.Infile[i],".asp"))){
      RootIndex = i-1;
      break;
    }
  }
  if(RootIndex == 0){
    printf("WARNING: input file does not follow .asp convention.\n");
    RootIndex = strlen(RunMode.Infile)-1;
  }

  strcpy(FitsFile,RunMode.Infile);

  if(fits_open_file(&Fin, FitsFile, READONLY, &status)){
    printf("Error opening FITS file %s !!!\n", FitsFile);
    exit(1);
  }

  printf("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("NEW VERSION!!!!!!!!!!!!!!!\n");
  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
  
  /* Read in values for header variables */
  if(ReadASPHdr(&InHdr, Fin) < 0){
    printf("%s> Unable to read Header from file %s.  Exiting...\n",ProgName,
	   RunMode.Infile);
    exit(2);
  }
  printf("\n==========================\n");
  printf("ASP FITS Header %s\n",InHdr.gen.HdrVer);
  printf("==========================\n\n");fflush(stdout);
  
  printf("Input file:  %s\n\n",RunMode.Infile);fflush(stdout);

  printf("PSR %s:\n",InHdr.target.PSRName);
  printf("--------------\n\n");
  printf("Centre Frequency: %6.1lf MHz\n\n",InHdr.obs.FSkyCent);fflush(stdout);

  /* Add up total number of dumps between all files */
  RunMode.NDumps = 0;

  /* Get number of HDUs in fits file */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  
  if(!strcmp(InHdr.gen.HdrVer,"Ver1.0")){
    RunMode.NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				   many non-data tables we will be using */
  }
  else if(!strcmp(InHdr.gen.HdrVer,"Ver1.0.1")){
    RunMode.NDumps = (NumHDU-3)/2;
  }
  else{
    printf("Do not recognize FITS file version number in header.\n");
    printf("This header %s. Exiting...\n",InHdr.gen.HdrVer);fflush(stdout);
    exit(3);
  }
  
 
  /* Now check last dump to make sure it wrote properly if scan 
     was ctrl-c'd -- if so, reduce NDumps by 1 to avoid disaster */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  fits_movabs_hdu(Fin,NumHDU,&hdutype, &status);
  /* find NPtsProf */
  fits_get_num_rows(Fin, &NPtsProf, &status);status=0;    
  

  //  if ((RunMode.NBins!=(int)NPtsProf)){
  if ((InHdr.redn.RNBinTimeDump!=(int)NPtsProf)){
    printf("Warning:  last dump did not write cleanly.\n");
    printf("Reducing number of dumps read in file by one...\n\n");
    fflush(stdout);
    RunMode.NDumps--;
  }
  
  /* Get and sort out command line options */
  if(GetOptions(&RunMode, &CalMode, Cmd, &InHdr) < 0){
    printf("Command line parsing failed.  Exiting...\n");fflush(stdout);
    exit(4);
  } 

  if(GetChans(&InHdr, Cmd, &RunMode) < 0) {
    printf("Could not parse channel options.  Exiting...\n");
    fflush(stdout);
    exit(5);
  }

  if(RunMode.Verbose){
    printf("AddChans = %d\n",Cmd->AddChansP);fflush(stdout);
    printf("NChan = %d\n",InHdr.obs.NChan);fflush(stdout);
    printf("NOutChans = %d\n",RunMode.NOutChans);fflush(stdout);
    printf("AddDumps = %d\n",RunMode.AddDumps);fflush(stdout);
    printf("NDumps = %d\n",RunMode.NDumps);fflush(stdout);
  }
  
  printf("Number of input channels:  %d\n",InHdr.obs.NChan);
  printf("Number of input dumps:     %d\n\n",RunMode.NDumps);

  printf("Number of output channels: %d\n",RunMode.NOutChans);
  printf("Number of output dumps:    %d\n\n",RunMode.NOutDumps);
  fflush(stdout);

  /* Set some output header variables */
  memcpy(&OutHdr, &InHdr, sizeof(struct ASPHdr));
  strcpy(OutHdr.gen.HdrVer,"Ver1.0.1");
  OutHdr.redn.RNTimeDumps = RunMode.NOutDumps;
  OutHdr.obs.NChan = RunMode.NOutChans;
  OutHdr.redn.RNBinTimeDump = RunMode.NBinsOut;
  for(i_chan_out=0;i_chan_out<OutHdr.obs.NChan;i_chan_out++) {
    OutChanIndex = (int)(((double)(RunMode.FirstChanAdd[i_chan_out] + 
				   RunMode.LastChanAdd[i_chan_out]))/2.);
    OutHdr.obs.ChanFreq[i_chan_out]=InHdr.obs.ChanFreq[OutChanIndex];
    printf("ChanfreqOut = %lf\n",OutHdr.obs.ChanFreq[i_chan_out]);
  }
  /* Set up output fits file stuff */
  NFilesOut = (int)(((float)RunMode.NOutDumps)/MAXDUMPS)+1;
  Fout = (fitsfile **)malloc(NFilesOut*sizeof(fitsfile));
  FileOutNo=-1;


  /* malloc pointers to input data arrays -- will use these to construct 
     Stokes parameters */
  ASquared    = (double **)malloc(InHdr.obs.NChan*sizeof(double));
  BSquared    = (double **)malloc(InHdr.obs.NChan*sizeof(double));
  ReAconjB    = (double **)malloc(InHdr.obs.NChan*sizeof(double));
  ImAconjB    = (double **)malloc(InHdr.obs.NChan*sizeof(double));
  SampleCount = (int    **)malloc(InHdr.obs.NChan*sizeof(int));

  /* malloc Output profiles */
  InputProfs = (struct StdProfs *)malloc(InHdr.obs.NChan*
					 sizeof(struct StdProfs));
  //  OutputProfs = (struct StdProfs *)malloc(InHdr.obs.NChan*
  OutputProfs = (struct StdProfs *)malloc(RunMode.NOutChans*
					  sizeof(struct StdProfs));  
  /*  AddChansProfs = (struct StdProfs *)malloc(RunMode.NumEffChans*
      sizeof(struct StdProfs)); */
  /* malloc subheaders (specific to each dump) */
  SubInHdr  = (struct SubHdr *)malloc(RunMode.NDumps*sizeof(struct SubHdr)); 


  /* Parse omission file and set up omission scans etc. */
  if (GetOmit(&InHdr, Cmd, &RunMode) < 0){
    printf("Unable to read or parse omit file %s. Exiting...\n",Cmd->Omitfile);
    fflush(stdout);
    exit(6);
  }

  /* malloc array of Cal Factors */
  for(i=0;i<NCHMAX;i++)
    JyPerCount[i] = (double *)malloc(4*sizeof(double));

  /* Get calibration factors (in Jy/count) from cal file, if given. */
  if(ReadCal(&InHdr, &RunMode, &CalMode, JyPerCount) < 0) {
    printf("Cannot perform calibration.  Exiting...\n");fflush(stdout);
    exit(7);
  }


  /* malloc header lines for ascii output */
  for(i_chan_in=0;i_chan_in<InHdr.obs.NChan;i_chan_in++){
   if( ((HeadLine[i_chan_in] = (char *)malloc(128)) == NULL) || 
	((AddChansHead[i_chan_in] = (char *)malloc(128)) == NULL) || 
	((OutputHead[i_chan_in] = (char *)malloc(128)) == NULL) ){
      printf("HeadLine malloc'ing failed for i_chan_in = %d\n",i_chan_in);
      fflush(stdout);
      exit(5);
    }
  }  

  /* Read in standard profile and rotate it to zero phase */
  if (ReadStd(&RunMode, &StdRunMode, &StdProfile, &TotWgt) < 0 ){
    printf("Could not read standard profile correctly.  Exiting...\n");
    fflush(stdout);
    exit(6);
  } 

  /* Read in ThetaBB values and Mueller matrix */
  if (RunMode.ThetaBBFlag){
    printf("Applying ThetaBB correction...\n");fflush(stdout);
    if(GetThetaBB(&RunMode, &InHdr) < 0){
      printf("Could not read in ThetaBB values correctly.  Exiting...\n");
      fflush(stdout);
      exit(7);
    }
  }
  if (Cmd->MuellerfileP){
    printf("Applying Mueller Matrix correction...\n");fflush(stdout);
    if(GetMueller(Cmd->Muellerfile, RunMode.MM, &InHdr) < 0){
      printf("Could not read in Mueller Matrix from file.  Exiting...\n");
      fflush(stdout);
      exit(8);
    }
  }

  /* Read in polyco.dat file and get polyco structure if requested on 
     command line */
  if(Cmd->PolyfileP){
    if(Cmd->PolyfileC == 0) sprintf(Cmd->Polyfile,"polyco.dat");
    printf("Polyco file name:  %s\n",Cmd->Polyfile);
    StartMJD = (double)InHdr.obs.IMJDStart + 
      ((double)InHdr.obs.StartTime)/86400.;
    /* malloc Polyco structure to number of dumps */
    Polycos = (struct Polyco *)malloc(MAX_PC_SETS*InHdr.obs.NChan*
				      sizeof(struct Polyco));    
    for(i_chan_in=0;i_chan_in<InHdr.obs.NChan;i_chan_in++){
      if((n_poly=GetPoly(Cmd->Polyfile, RunMode.Source, 
			 &Polycos[i_chan_in*MAX_PC_SETS], 
			 InHdr.obs.ChanFreq[i_chan_in], 
			 StartMJD)) < 1) {
	printf("Could not find polycos for all input profiles of \n");
	printf("PSR %s given as input.  Will not shift profiles.\n",
	       RunMode.Source);
	Cmd->PolyfileP=0;
      }
    }
    printf("Polycos successfully found.\n\n");
  }
  
  if(Cmd->PAngleOffP)
    printf("Will NOT fit out parallactic angle changes.\n");fflush(stdout);

  /* Move to the first data table HDU in the fits file */
  if(!strcmp(InHdr.gen.HdrVer,"Ver1.0"))
    fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT0", 0, &status);
  else if (!strcmp(InHdr.gen.HdrVer,"Ver1.0.1"))
    fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF0", 0, &status);

  /* Get the current HDU number */
  fits_get_hdu_num(Fin, &NFirstTable);

  if(RunMode.Verbose) {
    printf("\nPSR %s\n",RunMode.Source);
    printf("-----------\n");
    printf("DM = %f\n",InHdr.obs.DM);fflush(stdout);
  }

  
  for (i_dump_out=0; i_dump_out<RunMode.NOutDumps; i_dump_out++){
    /* If we are in the last Output dump, max dump to add is actual last 
       dump so we don't go over and add things that don't exist */

    MinAddDump = i_dump_out*RunMode.AddDumps; 
    if(i_dump_out == RunMode.NOutDumps-1){
      MaxAddDump = RunMode.NDumps; // not NDumps-1 to include last dump
    }
    else{
      MaxAddDump = MinAddDump + RunMode.AddDumps;
    }
        
  
    for (i_chan_out=0; i_chan_out<RunMode.NOutChans; i_chan_out++){
      /* Zero out each part of output profile, once for each output channel */
      FZero(OutputProfs[i_chan_out].rstds,NBINMAX); 
      FZero(OutputProfs[i_chan_out].rstdq,NBINMAX);
      FZero(OutputProfs[i_chan_out].rstdu,NBINMAX); 
      FZero(OutputProfs[i_chan_out].rstdv,NBINMAX); 
      FZero(OutputProfs[i_chan_out].stdlin,NBINMAX); 
      FZero(OutputProfs[i_chan_out].stdphi,NBINMAX); 
      FZero(OutputProfs[i_chan_out].stdphierr,NBINMAX);
    }
    

   /*** Loop over each integration ***/
      
    /* Now loop over beginning and end dumps going into each output dump */
    for (i_dump_in=MinAddDump; i_dump_in<MaxAddDump; i_dump_in++){
	
	
      /* move to next dump's data */
      if(!strcmp(InHdr.gen.HdrVer,"Ver1.0")){
	fits_movabs_hdu(Fin, NFirstTable+i_dump_in, &hdutype,&status);
	/* find NPtsProf */
	fits_get_num_rows(Fin, &NPtsProf, &status);status=0; 
      }
      else if(!strcmp(InHdr.gen.HdrVer,"Ver1.0.1")){
	/* if we've reached the end of the FITS file then increase FileNo */
	fits_movabs_hdu(Fin, NFirstTable+(i_dump_in)*2+1,&hdutype,
			&status);
	/* find NPtsProf */
	fits_get_num_rows(Fin, &NPtsProf, &status);status=0; 
	fits_movrel_hdu(Fin, -1, NULL, &status);
      }
	

      /* Exclude badly written dumps (usually due to ctrl-c)*/
      if (RunMode.NBins != (int)NPtsProf){
	if(RunMode.OldFits) { /* total hack */
	  printf("\nNumber of bins in header does not agree with header value.\n");
	  printf("\nSkipping rest of file...\n");fflush(stdout);
	  RunMode.NDumps = i_dump_in;
	  break; 
	}
	else{
	  printf("Warning: Number of bins in header (%d) does not agree with number \n",
		 RunMode.NBins);
	  printf("of bins in profiles (%d) for dump %d! Exiting...\n",
		 (int)NPtsProf, i_dump_in);fflush(stdout);
	  exit(8);
	}
      }
      /* Error check to make sure that Number of output 
	 bins < Number of input bins */
      if(RunMode.NBinsOut > (int)NPtsProf){
	printf("\n%s: Number of output bins (%d) must be smaller than\n",
	       ProgName,RunMode.NBinsOut);
	printf("number of input bins (%d) -- dump %d.  Skipping this dump...\n",
	       (int)NPtsProf, i_dump_in);fflush(stdout);
	exit(9);
      }
	
	
      /* Read in data arrays */
      ReadASPData(&InHdr, &SubInHdr[i_dump_in], &RunMode, Fin, i_dump_in,
		  NPtsProf, ASquared, BSquared, ReAconjB, ImAconjB, 
		  SampleCount, HeadLine);
	
	
      /* Now we have arrays for each channel for this dump */
      /* Create Stokes parameters for each profile in each channel */
	
      /*** Loop over each channel ***/
	
      /* Do one (output) channel at a time */
      for (i_chan_out=0; i_chan_out<RunMode.NOutChans; i_chan_out++){
            
      	/* Now loop over beginning and end input channels of each 
	   output channel */
	/* Note that we put "<=" here since we want the range to be 
	   inclusive */
	for (i_chan_in=RunMode.FirstChanAdd[i_chan_out]; 
	     i_chan_in<=RunMode.LastChanAdd[i_chan_out]; i_chan_in++){ 
	  
	  if(RunMode.Verbose)
	    printf("\n Channel number %d = %f MHz\n",i_chan_in,
		   InHdr.obs.ChanFreq[i_chan_in]);fflush(stdout);

	  /* Make Stokes profs, shift phases accordingly -- EVEN IF we are 
	     mitting the particular scan, since it may be used to designate 
	     the phase and time for the added profile; i.e. we don't want
	     the phase of the output added profile to be off by any amount.  */

	  /* Construct Stokes parameters */
	  MakeStokes(&InHdr, &RunMode, &InputProfs[i_chan_in], 
		     ASquared[i_chan_in], BSquared[i_chan_in],
		     ReAconjB[i_chan_in], ImAconjB[i_chan_in],
		     JyPerCount[CalMode.CalIndex[i_chan_in]]);
	  
	  /* Shift phases by apprpriate amounts if requested on command line 
	     -- do for each dump and channel */
	  if(Cmd->PolyfileP){
	    if(i_dump_in==0 && i_chan_in==0) {
	      printf("Applying polyco-based phase shifts...\n\n");
	      fflush(stdout);
	    }
	    if(PhaseShift(&Polycos[i_chan_in*MAX_PC_SETS],n_poly, 
			  &InputProfs[i_chan_in], &RunMode,
			  &InHdr, &SubInHdr[i_dump_in], i_chan_in) < 0) {
	      printf("Unable to shift profile phases.  Exiting...\n");
	      fflush(stdout);
	      exit(11); 
	    }
	  }
	  

	  /* Test to see if omitted.  If yes, can skip lots of stuff */
	  if(!RunMode.OmitFlag[i_dump_in*InHdr.obs.NChan + i_chan_in]) {

	    /* Create filename for un-added profile files */
	    sprintf(Stokesfile,"%s.%4.4d.%4.4d.prof.asc",RunMode.OutfileRoot,
		    (int)(InHdr.obs.ChanFreq[i_chan_in]),i_dump_in);
	    
	    if(RunMode.Verbose)
	      printf("Stokesfile  = %s\n",Stokesfile);fflush(stdout);
	    
	    
	    /* Bin Down(if input on command line) and Baseline subtract */
	    if(RunMode.BinDown)
	      BinDown(&RunMode, &InputProfs[i_chan_in], &StokesProfs);
	    else
	      memcpy(&StokesProfs, &InputProfs[i_chan_in], 
		     sizeof(struct StdProfs));
	    MakePol(&RunMode, RunMode.NBinsOut, &StokesProfs);
	    
	    /* Write Stokes parameters, linear polarization, 
	       and position angle to file if desired */
	    if (Cmd->WriteAllP)
	      WriteStokes(&RunMode, &StokesProfs, HeadLine[i_chan_in], 
			  Stokesfile);
	    
	    /* Correct for phase offset between polarizations */
	    if(Cmd->ThetaBBfileP)
	      FitThetaBB(&RunMode, &InHdr, &InputProfs[i_chan_in], 
			 i_chan_in, i_dump_in); 
	    
	    /* Fit out the Mueller Matrix */
	    if (Cmd->MuellerfileP)
	      FitMueller(&RunMode, &InHdr, &InputProfs[i_chan_in], i_chan_in);
	    
	    /* Correct for parallactic angle -- an option!! */
	    if (!Cmd->PAngleOffP)
	      FitAngle(&RunMode, &InHdr, &SubInHdr[i_dump_in], 
		       &InputProfs[i_chan_in]);  
	    
	    
	    /* Update highest-SNR profile to use as standard prof if 
	       none provided */
	    InputProfs[i_chan_in].SNR = StokesProfs.SNR;
	    
	    /*** print RMS values of input profiles to file for checks ***/
	    
	    
	    /* Cumulatively add Output profiles */

	    if (RunMode.Verbose)
	      printf("\nWill now add channels together...\n");fflush(stdout);
	    
	    
	    for (i_bin=0; i_bin<RunMode.NBins; i_bin++){
	      /**** add maybe routine to do simple cumulative addition ****/
	      OutputProfs[i_chan_out].rstds[i_bin] += 
		InputProfs[i_chan_in].rstds[i_bin];
	      OutputProfs[i_chan_out].rstdq[i_bin] += 
		InputProfs[i_chan_in].rstdq[i_bin];
	      OutputProfs[i_chan_out].rstdu[i_bin] += 
		InputProfs[i_chan_in].rstdu[i_bin];
	      OutputProfs[i_chan_out].rstdv[i_bin] += 
		InputProfs[i_chan_in].rstdv[i_bin];
	    }  /**** end add ***/

	  } /* end of if (OmitFlag) */
	  
	
	} /* end of i_chan_in loop */
	
      
       
     
      
      }  /* end i_chan_out loop */
      
    } /* end i_dump_in loop */
    

    /** Set Middle time stamp in SubOutHdr for this output dump **/
    SubOutDumpIndex = (int)(floor((double)(MinAddDump + MaxAddDump - 1))/2.);
    SubOutHdr.DumpMiddleSecs = SubInHdr[SubOutDumpIndex].DumpMiddleSecs;

    /** Now set SubHdr entries for phase and period **/
    /** It is the central dump, and central channel **/
    /** Need to do it here because need to have had read all dumps from 
	this set in order to have data in hand **/
    for (i_chan_out=0; i_chan_out<RunMode.NOutChans; i_chan_out++){
      for (i_chan_in=RunMode.FirstChanAdd[i_chan_out]; 
	   i_chan_in<=RunMode.LastChanAdd[i_chan_out]; i_chan_in++){ 
	SubOutChanIndex = (int)(((double)(RunMode.FirstChanAdd[i_chan_out] + 
					  RunMode.LastChanAdd[i_chan_out]))/2.);
	SubOutHdr.DumpRefPhase[i_chan_out] = 
	  SubInHdr[SubOutDumpIndex].DumpRefPhase[SubOutChanIndex];
	SubOutHdr.DumpRefPeriod[i_chan_out] = 
	  SubInHdr[SubOutDumpIndex].DumpRefPeriod[SubOutChanIndex];
      }
     /* Divide profiles by proper denominator for averaging profiles */
      i_out=i_dump_out*RunMode.NOutChans+i_chan_out;
      denom = RunMode.TotScans[i_out]-RunMode.TotOmit[i_out];
      /* if for some reason # scan < # omissions (shouldn't happen), 
	 then barf */
      if(denom<0) {
	printf("Not sure why, but more omissions than scans went into ");
	printf("this profile.  Let Rob know.  Exiting...\n");
	fflush(stdout);
	exit(12);
      }
      /* If all input scans going into given output profile are omitted, set
	 entire profile to zero for now */
      if(denom==0){ 
	FZero(OutputProfs[i_chan_out].rstds,NBINMAX); 
	FZero(OutputProfs[i_chan_out].rstdq,NBINMAX);
	FZero(OutputProfs[i_chan_out].rstdu,NBINMAX); 
	FZero(OutputProfs[i_chan_out].rstdv,NBINMAX); 
      }
      /* Otherwise, average profile! */
      else {
	for (i_bin=0; i_bin<RunMode.NBins; i_bin++){
	  OutputProfs[i_chan_out].rstds[i_bin] /= (float)denom;
	  OutputProfs[i_chan_out].rstdq[i_bin] /= (float)denom;
	  OutputProfs[i_chan_out].rstdu[i_bin] /= (float)denom;
	  OutputProfs[i_chan_out].rstdv[i_bin] /= (float)denom;
	}
      }
      /****** bin down and make pol'n ******/
      if(RunMode.BinDown){
	memcpy(&TempOutProfs, &OutputProfs[i_chan_out], 
	       sizeof(struct StdProfs));
	BinDown(&RunMode, &TempOutProfs, &OutputProfs[i_chan_out]);
      }
      MakePol(&RunMode, RunMode.NBinsOut, &OutputProfs[i_chan_out]);
    }
    

    /* BEGIN WRITING TO FITS FILE */
    
    if(RunMode.Verbose)
      printf("WRITING OUTPUT FITS FILE\n");fflush(stdout);
    /* test to see what new Dump information is */
    if(RunMode.Verbose){
      printf("\n\nOutput Dump %d:  TIME OF DUMP = %lf\n",i_dump_out,
	     SubOutHdr.DumpMiddleSecs);
      printf("          CHANNEL (MHz)   REF. PHASE   REF. PERIOD (s)\n");
      printf("          -------------   ----------   ---------------\n");
      for(i_chan_out=0;i_chan_out<OutHdr.obs.NChan;i_chan_out++) 
	printf("          %13.1lf%13.8lf%18.11lf\n", 
	       OutHdr.obs.ChanFreq[i_chan_out],
	       SubOutHdr.DumpRefPhase[i_chan_out],
	       SubOutHdr.DumpRefPeriod[i_chan_out]);
      printf("\n\n");fflush(stdout);
    }
    
    if(i_dump_out%MAXDUMPS == 0) {
      FileOutNo++;
      
      strcpy(FitsFileOut,"\0");
      if(NFilesOut > 1){
	sprintf(FitsFileOut,"%s.%d.stokes.fits",RunMode.OutfileRoot,FileOutNo);
      }
      else{
	sprintf(FitsFileOut,"%s.stokes.fits",RunMode.OutfileRoot);
      }
      
      printf("Writing file %s\n\n",FitsFileOut);fflush(stdout);
      
      if(fits_create_file(&Fout[FileOutNo], FitsFileOut, &fitsstatus)) {
	printf("Error opening FITS file %s!!!\n",
	       FitsFileOut);fflush(stdout);
	exit(10);
      }
      
      if(WrtASPHdr(&OutHdr, Fout[FileOutNo]) < 0) {
	printf("Error writing ASP header structure!\n");fflush(stdout);
	exit(11);
      }
    }    
    
    /* only include nonzero (i.e. not all-omitted) added scans */ 
    /*    if (denom == 0 && Cmd->NoBadP) { 
      printf("Output scan %d is excluded from final output.\n",
	     i_dump_out);
    }
    else{  */
    if(WrtASPStokes(OutHdr, SubOutHdr, Fout[FileOutNo], i_dump_out, 
		    OutputProfs, Cmd->OmitfileP, &RunMode) < 0) {
      printf("Cannot write data tables to fits file. Exiting...\n");
      exit(12);
    }
      /*  }  */
      
  } /* end of i_dump_out loop*/
  
  free(OutputProfs);
  free(InputProfs);
  //  free(TempOutProfs);
  
  
  /* Close FITS files */
  fits_close_file(Fin, &fitsstatus);
  for(i=0;i<NFilesOut;i++)
    fits_close_file(Fout[i], &fitsstatus);

  printf("Completed successfully.\n\n");fflush(stdout);

  PrintLog(&RunMode, &InHdr, Cmd);

  exit(0);
}



void PrintLog(struct RunVars *RunMode, struct ASPHdr *Hdr, Cmdline *Cmd) 
{

  int    i_chan_in, i_chan_out, i_dump_in, i_dump_out;
  int    OutChanIndex;
  double ChanFreq;

  FILE   *Fcheck;
  /*  sprintf(testfile,"test_denom.dat");
  Fcheck=fopen(testfile, "a");
  fprintf(Fcheck, "ADDCHANS DENOMINATOR -- DUMP %d:\n\n",nscan);
  fprintf(Fcheck,"Combined Channel %d: \n", i);
  fprintf(Fcheck," %d      --->  %d\n",
	  RunMode->MinChans2Add[i],RunMode->MaxChans2Add[i]);
  fprintf(Fcheck," %6.1lf  --->  %6.1lf\n", 
	  Hdr->obs.ChanFreq[RunMode->MinChans2Add[i]], 
	  Hdr->obs.ChanFreq[RunMode->MaxChans2Add[i]]);
  fprintf(Fcheck, "NumChans2Add = %d\n",RunMode->NumChans2Add[i]);
  fprintf(Fcheck, "zapchans = %d\n",zapchans);
  fprintf(Fcheck, "omitcount = %d\n",omitcount);
  fprintf(Fcheck, "DENOMINATOR = %d\n\n\n",RunMode->NumChans2Add[i] - zapchans - omitcount);
  fclose(Fcheck);  */
  
  

  Fcheck=fopen("check_chans.dat","a");
  fprintf(Fcheck,"Original Channels: \n\n");
  for(i_chan_in=0; i_chan_in<Hdr->obs.NChan; i_chan_in++){
    fprintf(Fcheck,"  %6.1lf  ",Hdr->obs.ChanFreq[i_chan_in]);
  }
  fprintf(Fcheck,"\n\n\n");
  fprintf(Fcheck,"Combined Channels (Total %d):\n\n", RunMode->NOutChans);
  for(i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
    OutChanIndex = (int)(((double)(RunMode->FirstChanAdd[i_chan_out] + 
				   RunMode->LastChanAdd[i_chan_out]))/2.);
    ChanFreq = Hdr->obs.ChanFreq[OutChanIndex];
    fprintf(Fcheck,"Output Channel %d -- %6.1lf MHz: \n", i_chan_out, ChanFreq);
    fprintf(Fcheck," %d      --->  %d\n",
	    RunMode->FirstChanAdd[i_chan_out],
	    RunMode->LastChanAdd[i_chan_out]);
    fprintf(Fcheck," %6.1lf  --->  %6.1lf\n\n", 
	    Hdr->obs.ChanFreq[RunMode->FirstChanAdd[i_chan_out]], 
	    Hdr->obs.ChanFreq[RunMode->LastChanAdd[i_chan_out]]);
    /*      fprintf(Fcheck,"Total Scans: %d\n",RunMode->TotScans[nscan*RunMode->NOutChans + i_chan_out]);
	    fprintf(Fcheck,"Total Omit:  %d\n\n",RunMode->TotOmit[nscan*RunMode->NOutChans + i_chan_out]);  */
  }
  fclose(Fcheck);


  if(Cmd->OmitfileP){
      
    Fcheck=fopen("check_omit.dat", "w");
    fprintf(Fcheck, "  Omit matrix for each Input Dump/Channel combination (omit==1):\n\n");
    fprintf(Fcheck, "  ");
    for (i_chan_in=0; i_chan_in<Hdr->obs.NChan; i_chan_in++)
      fprintf(Fcheck, " %6.1lf ",Hdr->obs.ChanFreq[i_chan_in]);
    fprintf(Fcheck, "\n\n");
      
    for (i_dump_in=0; i_dump_in<RunMode->NDumps; i_dump_in++){
      fprintf(Fcheck, "%d  ",i_dump_in);
      for (i_chan_in=0; i_chan_in<Hdr->obs.NChan; i_chan_in++){
	fprintf(Fcheck,"   %d    ",RunMode->OmitFlag[i_dump_in*Hdr->obs.NChan + i_chan_in]);
      }
      fprintf(Fcheck, "\n");
    }
        fprintf(Fcheck, "\n\n\n");
    
    //    Fcheck=fopen("test_omit_totals.dat", "w");
    fprintf(Fcheck, "  Total Input Scans in each Output Dump/Channel combination:\n\n");
    fprintf(Fcheck, "   ");
      
      
    for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
      /* OutChan is average channel */
      OutChanIndex = (int)(((double)(RunMode->FirstChanAdd[i_chan_out] + 
				     RunMode->LastChanAdd[i_chan_out]))/2.);
      fprintf(Fcheck, " %6.1lf ",Hdr->obs.ChanFreq[OutChanIndex]);
    }
    fprintf(Fcheck, "\n\n");
      
    for (i_dump_out=0; i_dump_out<RunMode->NOutDumps; i_dump_out++){
      fprintf(Fcheck, "%d  ",i_dump_out);
      for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
	fprintf(Fcheck,"   %d    ",
		RunMode->TotScans[i_dump_out*RunMode->NOutChans + i_chan_out]);
      }
      fprintf(Fcheck, "\n");
    }
    fprintf(Fcheck, "\n\n\n");
    fprintf(Fcheck, "  Total Omissions from each Output Dump/Channel combination:\n\n");
    for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
      /* OutChan is average channel */
      OutChanIndex = (int)(((double)(RunMode->FirstChanAdd[i_chan_out] + 
				     RunMode->LastChanAdd[i_chan_out]))/2.);
      fprintf(Fcheck, " %6.1lf ",Hdr->obs.ChanFreq[OutChanIndex]);
    }
    fprintf(Fcheck, "\n\n");
      
    for (i_dump_out=0; i_dump_out<RunMode->NOutDumps; i_dump_out++){
      fprintf(Fcheck, "%d  ",i_dump_out);
      for (i_chan_out=0; i_chan_out<RunMode->NOutChans; i_chan_out++){
	fprintf(Fcheck,"   %d    ",
		RunMode->TotOmit[i_dump_out*RunMode->NOutChans + i_chan_out]);
      }
      fprintf(Fcheck, "\n");
    }
    fclose(Fcheck);

  }
	
  return;
	
}


//Print
