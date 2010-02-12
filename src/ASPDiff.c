/* Takes output stokes profile fits file and subtracts each profile therein 
   from an input ASCII profile, after rotating and scaling 

   * --R. Ferdman,  10 November 2009 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "DiffCmdLine.h" 
#include "ASPCommon.h"

int main(int argc, char *argv[])
{
  

  int NDumps, bin[NBINMAX], NAscBins;
  int NFirstTable, NumHDU, hdutype;
  struct ASPHdr Hdr;
  struct SubHdr StokesSubHdr;

  struct StdProfs AscProfile, *Profile, *DiffProfile;
  int OutRootIndex=0, LastSlashIndex=0, TotalRootIndex=0;
  char FitsFile[128], OutRoot[128], ProgName[32];
  char DiffFile[64];
  char AscHeader[256];
  fitsfile *Fstokes;
  fitsfile *Fdiff;

  int ngood;
  float  profs[NBINMAX],amps[NBINMAX], phas[NBINMAX];
  float  Shift,EShift,SNR,ESNR,b,errb;
  double ByAngle=0.;

  int i, i_bin, nscan, status=0;
  int fitsstatus=0;
  long NPtsProf=0;
  double RA, Dec;
  char NObs[3];

  double SBase,Srms,VBase,Vrms,UBase,QBase,Urms,Qrms;
  double Duty,SPeak;
  double FinalMask[NBINMAX];
  int spk;
  Cmdline         *Cmd;

  struct RunVars  RunMode;

                                                                             
 /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /* Store program name */
  strcpy(ProgName, argv[0]);


  /* read in ascii profile */
  if (Cmd->AscFileP){
    
    if (ReadASPAsc(Cmd->AscFile, AscHeader, bin, &AscProfile, &NAscBins) < 0) {
      printf("Error in reading file %s.\n",Cmd->AscFile);
      fflush(stdout);
      exit(1);
    }
    // strcpy(RunMode.Source,Cmd->Pulsar);
    RemoveBase(&RunMode, NAscBins, &AscProfile);

    /* cprofc the ascii profile */
    cprofc(AscProfile.rstds,NAscBins, AscProfile.stdamps, AscProfile.stdphas);
    
    
    memcpy(amps,AscProfile.stdamps,sizeof(float)*NBINMAX);
    memcpy(phas,AscProfile.stdphas,sizeof(float)*NBINMAX);
    

  }
  else{
    printf("Must supply an ASCII profile from which to subtract data ");
    printf("profiles.  Exiting...\n");
    exit(2);
  }



  strcpy(FitsFile, Cmd->Infile);
  
  LastSlashIndex = -1;
  for(i=0;i<strlen(FitsFile);i++){
    if(!strncmp(&FitsFile[i],"/",1))
      LastSlashIndex = i;
  }
  /* Get output file root name from input file name */
  for (i=0;i<strlen(FitsFile);i++){
    if(!(strcmp(&FitsFile[i],"stokes.fits"))){
      OutRootIndex = i-1;
      break;
    }
  }

  if(OutRootIndex == 0){
    printf("WARNING: input file name does not follow .stokes.fits convention.\n");
    sprintf(OutRoot,"ASPOut");
  }
  else {
    TotalRootIndex = OutRootIndex-LastSlashIndex-1;
    strcpy(OutRoot,"\0");
    strncpy(OutRoot,&FitsFile[LastSlashIndex+1],TotalRootIndex);
    strcpy(&OutRoot[TotalRootIndex],"\0");
  }

  /* Malloc fitsfile descriptor to number of data files */
  //  Fstokes  = (fitsfile **)malloc(Cmd->NFiles*sizeof(fitsfile));

  /* Open fits data file */

  if(fits_open_file(&Fstokes, FitsFile, READONLY, &status)){
    printf("Error opening FITS file %s !!!\n", FitsFile);
    exit(1);
  }
  
  /* Read in values for header variables */
  if(ReadASPHdr(&Hdr, Fstokes) < 0){
    printf("%s> Unable to read Header from file %s.  Exiting...\n",
	   ProgName,FitsFile);
    exit(2);
  }

  printf("\n==========================\n");
  printf("ASP FITS Header %s\n",Hdr.gen.HdrVer);
  printf("==========================\n\n");fflush(stdout);
  
  printf("Input file:  %s\n\n",FitsFile);fflush(stdout);

  printf("PSR %s:\n",Hdr.target.PSRName);
  printf("--------------\n\n");
  printf("Centre Frequency: %6.1lf MHz\n\n",Hdr.obs.FSkyCent);fflush(stdout);
  //  if(Cmd->NoBaseP)
  //  printf("Baseline subtraction turned OFF.\n\n");fflush(stdout);

  
  NDumps = 0;
  /* Get number of HDUs in fits file */
  fits_get_num_hdus(Fstokes, &NumHDU, &status);
  
  /* Temporary -- need to write to and read this from Header */
  if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
    NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				 many non-data tables we will be using */
  }
  else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
    NDumps = (NumHDU-3)/2;
  }
  else{
    printf("Do not recognize FITS file version number in header.\n");
    printf("This header %s. Exiting...\n",Hdr.gen.HdrVer);fflush(stdout);
    exit(3);
  }
  //  if(Cmd->NFiles > 1)
  //    printf("Total number of dumps: %d\n",NDumps);
  
  printf("Number of channels:  %d\n",Hdr.obs.NChan);
  printf("Number of dumps:     %d\n\n",NDumps);

  Profile = (struct StdProfs *)malloc(Hdr.obs.NChan*sizeof(struct StdProfs));
  DiffProfile =(struct StdProfs *)malloc(Hdr.obs.NChan*sizeof(struct StdProfs));

  /* Conver Ra and Dec, mainly for use with parallactic angle calculation if needed */
  if(Hdr.target.RA > 0. && fabs(Hdr.target.Dec) > 0. ){
    RA = (Hdr.target.RA/24.)*TWOPI;
    Dec = Hdr.target.Dec*TWOPI/360.;
  }
  else{
    printf("\nRA and Dec are not in file.  Parallactic angle calculations\n");
    printf(" will be incorrect...\n");fflush(stdout);
  }
  /* Convert NObs to integer for possible use in calculating parallactic angle */
  sscanf(Hdr.obs.ObsvtyCode,"%s",NObs);

  /* Move to the first data table HDU in the fits file */
  if(!strcmp(Hdr.gen.HdrVer,"Ver1.0"))
    fits_movnam_hdu(Fstokes, BINARY_TBL, "STOKES0", 0, &status);
  else if (!strcmp(Hdr.gen.HdrVer,"Ver1.0.1"))
    fits_movnam_hdu(Fstokes, ASCII_TBL, "DUMPREF0", 0, &status);

  /* Get the current HDU number */
  fits_get_hdu_num(Fstokes, &NFirstTable);


  /* open new fits file for writing */
  sprintf(DiffFile,"%s.diff.fits", OutRoot);

  if(fits_create_file(&Fdiff, DiffFile, &fitsstatus)) {
    printf("Error opening FITS file %s!!!\n",DiffFile);fflush(stdout);
    exit(4);
  }

  if(Cmd->VerboseP){
    printf("Rotating profile(s) by %6.3lf radians or %6.2lf degrees...\n\n",
	 ByAngle, ByAngle*360./TWOPI);
  }

  printf("Writing differenced file %s\n\n",DiffFile);fflush(stdout);

  /* Write header info into file */
  if(WrtASPHdr(&Hdr, Fdiff) < 0) {
    printf("Error writing ASP header structure!\n");fflush(stdout);
    exit(5);
  }
  
  /* Set up RunMode structure to be compatible with WrtASPStokes */
  RunMode.AddChans = 0;
  RunMode.AddDumps = 0;
  RunMode.NScanOmit = 0;
  IZero(RunMode.DumpOmit,MAXOMIT);
  IZero(RunMode.ChanOmit,MAXOMIT);
  IZero(RunMode.ZapChan,Hdr.obs.NChan);
  RunMode.NBins = RunMode.NBinsOut = Hdr.redn.RNBinTimeDump;
  
  if(NAscBins != RunMode.NBins) {
    printf("\nBoth ASCII and FITS profiles MUST have same number of bins! ");
    printf("(%d != %d) Exiting...\n",RunMode.NBins,NAscBins);
    fflush(stdout);
    exit(3);
  }
  
  /* now start running through each dump */
  for(nscan=0;nscan<NDumps;nscan++){

    /* move to next dump's data */
    if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
      fits_movabs_hdu(Fstokes, NFirstTable+nscan, &hdutype, &status); 
    }
    else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
      fits_movabs_hdu(Fstokes,NFirstTable+(nscan)*2+1,&hdutype,&status);
      fits_get_num_rows(Fstokes, &NPtsProf, &status);status=0; 
      fits_movrel_hdu(Fstokes, -1, NULL, &status);
    }


    /* find NPtsProf */
    ReadASPStokes(&Hdr, &StokesSubHdr, Fstokes, NPtsProf, 
		  Profile, nscan, Cmd->VerboseP);


    for(i=0;i<Hdr.obs.NChan;i++){

    /*********** READING DONE.  BEGIN DIFFERENCE ************/
      
      /* First rotate and scale current fits profile */
      memcpy(profs,Profile[i].rstds,sizeof(float)*NBINMAX);

      /* Now fftfit to find shift required in second profile */
      fftfit_(profs,&amps[1],&phas[1],
	      &RunMode.NBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);  

      /* Now shift second profile -- convert Shift from bins to radians */
      
      if(!Cmd->NoRotP)
	ByAngle = (double)(-Shift/RunMode.NBins*TWOPI);
	RotateProf(&RunMode, &Profile[i], ByAngle);
      if(!Cmd->NoScaleP){
	for(i_bin=0;i_bin<RunMode.NBins;i_bin++){
	  Profile[i].rstds[i_bin] /= b;
	  Profile[i].rstdq[i_bin] /= b;
	  Profile[i].rstdu[i_bin] /= b;
	  Profile[i].rstdv[i_bin] /= b;
	} 
      }
      /* Now take the difference */
 
      MakePol(&RunMode,RunMode.NBins, &Profile[i]);
            
      /*********** ROTATION COMPLETED.  BEGIN WRITE TO FITS FILE ***********/
      

      Duty = DutyLookup(Hdr.target.PSRName);
      BMask(Profile[i].rstds,&Hdr.redn.RNBinTimeDump,&Duty,FinalMask);
      Baseline(Profile[i].rstds,FinalMask,&Hdr.redn.RNBinTimeDump,
	       &SBase,&Srms);
      Baseline(Profile[i].rstdv,FinalMask,&Hdr.redn.RNBinTimeDump,
	       &VBase,&Vrms);
      Baseline(Profile[i].rstdq,FinalMask,&Hdr.redn.RNBinTimeDump,
	       &QBase,&Qrms);
      Baseline(Profile[i].rstdu,FinalMask,&Hdr.redn.RNBinTimeDump,
	       &UBase,&Urms);
 

      SPeak =  FindPeak(Profile[i].rstds,&Hdr.redn.RNBinTimeDump,&spk);
      Profile[i].SNR = SPeak*Srms;
	
      if(!ArrayZero(Profile[i].rstds, RunMode.NBins)) { // i.e. good data
	for(i_bin=0;i_bin<RunMode.NBins;i_bin++){
	  DiffProfile[i].rstds[i_bin] = 
	    AscProfile.rstds[i_bin] - Profile[i].rstds[i_bin];
	  DiffProfile[i].rstdq[i_bin] = 
	    AscProfile.rstdq[i_bin] - Profile[i].rstdq[i_bin];
	  DiffProfile[i].rstdu[i_bin] = 
	    AscProfile.rstdu[i_bin] - Profile[i].rstdu[i_bin];
	  DiffProfile[i].rstdv[i_bin] = 
	    AscProfile.rstdv[i_bin] - Profile[i].rstdv[i_bin];
	} 
      }
      else{ // just replicate profile that is bad, i.e. zeroed-out
	for(i_bin=0;i_bin<RunMode.NBins;i_bin++){
	  DiffProfile[i].rstds[i_bin] = Profile[i].rstds[i_bin];
	  DiffProfile[i].rstdq[i_bin] = Profile[i].rstdq[i_bin];
	  DiffProfile[i].rstdu[i_bin] = Profile[i].rstdu[i_bin];
	  DiffProfile[i].rstdv[i_bin] = Profile[i].rstdv[i_bin];
	} 	
      }

    }

   if(WrtASPStokes(Hdr, StokesSubHdr, Fdiff, nscan, 
		    DiffProfile, 0, &RunMode) < 0) {
      printf("Cannot write data tables to fits file. Exiting...\n");fflush(stdout);
      exit(6);
    }
    
   if(Cmd->VerboseP) printf("\n");fflush(stdout);
  }


  //  for (i=0;i<Cmd->NFiles;i++)
  fits_close_file(Fstokes, &fitsstatus);
  fits_close_file(Fdiff, &fitsstatus);

  printf("Completed successfully.\n\n");fflush(stdout);

  exit(0);

}
