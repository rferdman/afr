/* Takes output stokes profile fits file and calculate TOAs using a 
 * user-supplied standard profile.  Based on pspmtoa routine.
 * 
 * --R. Ferdman,  18 July 2008 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "ToaCmdLine.h"  // change to a new Cmd file like ToaCmd.h
#include "ASPCommon.h"

int main(int argc, char *argv[])
{

  int i_bin, i_file, i_scan, i_chan, i_char; // counters
  int NFirstTable, NumHDU, hdutype, status=0; // FITS int's
  int RootIndex=0, LastSlashIndex=0; // filename indices
  int bin[NBINMAX], StdBins, NDumps=0; // profile info
  int NToa=0, TotToa=0, NWrtToa=0, TotWrtToa=0; // TOA counters: All TOAs, and Written-to-file TOAs (which may be different due to comman-line restrictions)
  int NOmit=0; // keep track of number of scans omitted from each file
  int ngood; // for fftfit
  int MJD0; // integer MJD for output TOA
  int freqwrite=0, mjdwrite=0, errwrite=0; // flags used to signify allowed freq and mjd ranges for TOAs

  long NPtsProf=0; // needs to ba a long for fits routines

  float ProfSum;
  float profs[NBINMAX],stdamps[NBINMAX], stdphas[NBINMAX], pha1; // for fftfit
  float Shift,EShift,SNR,ESNR,b,errb; // more fftfit variables

  double PhaseShift, AngleShift, PhaseChange; // calculated phase differences in two units
  double FracMJD, TOAErr; // Fractional MJD after shift, TOA error
  double RA, Dec; // not needed I don't think

  char FitsFile[128], Infile[64], Toafile[128]; // data file names
  char Headerline[256]; // header line in standard profile
  char ProgName[32]; // sotres this program name (ASPToa)
  char Source[8];
  char NObs[3]; // stores observatory code
  char TOA[32], tempMJD[9];
  char frac[24];
  char AllT2Flags[128];

  struct ASPHdr Hdr; // ASP header
  struct SubHdr StokesSubHdr; // SubHeader from *.stokes.fits files
  struct StdProfs StdProfile, *Profiles; // stores profile data

  FILE *Ftoa, *Fcheck; // output TOA file

  fitsfile *Fstokes; // file pointer for ASP *.stokes.fits file

  Cmdline *Cmd; // stores command line arguments
                              

                                               
  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /* Store program name */
  strcpy(ProgName, argv[0]);  

  /* Create an output file to check omissions if in vebose mode */
  if (Cmd->VerboseP || Cmd->CheckOmitP){
    if((Fcheck = fopen("check_omit.dat","w")) == 0)
      { printf("Cannot open check_omit.dat. Exiting...\n"); exit(1); }   
  }

  /** read in standard profile (ascii) and get phase info to 
      compare with each fits file scan **/

  if ( ReadASPAsc(Cmd->Template, Headerline, bin,  
		  &StdProfile, &StdBins) < 0) {
    printf("Error in reading file %s.\n",Cmd->Template);
    fflush(stdout);
    exit(1);
  }
  
  /* Now have standrad profile read in.  cprofc it: */
  cprofc(StdProfile.rstds,StdBins,StdProfile.stdamps,StdProfile.stdphas);

  /* Get amplitude and phase information from the standard profile */
  memcpy(stdamps,StdProfile.stdamps,sizeof(float)*NBINMAX);
  memcpy(stdphas,StdProfile.stdphas,sizeof(float)*NBINMAX);

  /* Rotate standard profile to zero phase unless user requests not to do so */
  pha1 = stdphas[1];
  if (!Cmd->NoZeroP){
    for(i_bin=1;i_bin<(StdBins/2)+1;i_bin++) 
      stdphas[i_bin] = fmod(stdphas[i_bin] -i_bin*pha1,TWOPI);
  }

  /** read it in so it only reads the total power and ignores lines with # 
      Maybe enforce 1st column=bin and second=total power, and further 
      columns are ignored **/

  /* Open output TOA file in order to write out TOAs line by line */

  if (Cmd->ToafileP) 
    strcpy(Toafile, Cmd->Toafile);
  else
    strcpy(Toafile, "toa.out");

  if((Ftoa = fopen(Toafile,"w"))==NULL){
    printf("Could not read file %s. Exiting...\n",Toafile);
    fflush(stdout);
    return -1;;
  }
  

  printf("\n=======\n");
  printf("ASPToa  %s\n",Hdr.gen.HdrVer);
  printf("=======\n\n");fflush(stdout);
     
  printf("Standard profile: %s\n\n",Cmd->Template);
  printf("Finding TOAs for %d data files.\n\n",Cmd->InfilesC);

  if (Cmd->NoZeroP){
    printf("Have chosen not to rotate standard profile to zero phase.\n");
  }
  else{
    printf("Rotated standard profile by %f radians to put it at ", pha1);
    printf("zero phase.\n");
  }

  if (Cmd->FreqRangeP) {
    if (Cmd->FreqRange[0] > Cmd->FreqRange[1]) {
      printf("Error:  first argument to -freq must be less than second ");
      printf("argument.  Exiting...");
      exit(2);
    }
    printf("Including only frequencies from %.1lf to %.1lf MHz\n",
	   Cmd->FreqRange[0], Cmd->FreqRange[1]);
  }
  if (Cmd->MJDRangeP) {
    if (Cmd->MJDRange[0] > Cmd->MJDRange[1]) {
      printf("Error:  first argument to -mjd must be less than second ");
      printf("argument.  Exiting...");
      exit(2);
    }
    printf("Including only TOAs with MJDs from %.1lf to %.1lf\n",
	   Cmd->MJDRange[0], Cmd->MJDRange[1]);
  }
  printf("\n+=========================================================================+\n\n");


  if(Cmd->Tempo2P){
    printf("Have chosen to output TOAs in tempo2 format\n\n");
    if(Cmd->MJDFlagP || Cmd->BEFlagP || Cmd->T2FlagsP){
      printf("Adding the following flags to the TOA line:\n");
      if(Cmd->MJDFlagP)
	printf("   -mjd <(shortened) MJD of each TOA>\n");
      if(Cmd->BEFlagP)
	printf("   -be <backend>\n");
      if(Cmd->T2FlagsP)
	printf("   %s\n",Cmd->T2Flags);
      printf("\n");
    }
  }
  else{

  }
  
  /**** FOR LOOP OVER N_FILES OR SOMETHING ****/

  for(i_file=0;i_file<Cmd->InfilesC;i_file++){ // over all input fits files

    NToa=0;  /* Initialize NToa for this file */
    NWrtToa=0;
    NOmit=0;

    strcpy(FitsFile, Cmd->Infiles[i_file]);

    /* Write file name in verbose-mode omit check file */
    if(Cmd->VerboseP || Cmd->CheckOmitP) 
      fprintf(Fcheck, "\n%s\n",Cmd->Infiles[i_file]);
  
    LastSlashIndex = -1;
    for(i_char=0;i_char<strlen(FitsFile);i_char++){
      if(!strncmp(&FitsFile[i_char],"/",1))
	LastSlashIndex = i_char;
    }
    /* Get output file root name from input file name */
    for (i_char=0;i_char<strlen(FitsFile);i_char++){
      if(!(strcmp(&FitsFile[i_char],"stokes.fits"))){
	RootIndex = i_char-1;
	break;
      }
    }
    strcpy(Infile,&FitsFile[LastSlashIndex+1]);

    if(RootIndex == 0){
      printf("WARNING: input file name does not follow .stokes.fits convention.\n");
    }

     /* Open fits data file, one at a time */

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

  
    printf("Input file:  %s\n",Infile);fflush(stdout);
     //   printf("Input file:  %s\n\n",FitsFile);fflush(stdout);

     printf("PSR %s\n",Hdr.target.PSRName);
     printf("Centre Frequency: %6.1lf MHz\n",Hdr.obs.FSkyCent);fflush(stdout);
  
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
  
    if(Cmd->VerboseP){
      printf("Number of channels:  %d\n",Hdr.obs.NChan);
      printf("Number of dumps:     %d\n\n",NDumps);
    }
    
    Profiles = (struct StdProfs *)malloc(Hdr.obs.NChan*sizeof(struct StdProfs));

    /* Convert Ra and Dec, mainly for use with parallactic angle 
       calculation if needed */
    if(Hdr.target.RA > 0. && fabs(Hdr.target.Dec) > 0. ){
      RA = (Hdr.target.RA/24.)*TWOPI;
      Dec = Hdr.target.Dec*TWOPI/360.;
    }
    else{
      printf("\nRA and Dec are not in file.  Parallactic angle calculations\n");
      printf(" will be incorrect...\n");fflush(stdout);
    }

    /* read in Observatoy code */
    sscanf(Hdr.obs.ObsvtyCode,"%s",NObs);

    /* Move to the first data table HDU in the fits file */
    if(!strcmp(Hdr.gen.HdrVer,"Ver1.0"))
      fits_movnam_hdu(Fstokes, BINARY_TBL, "STOKES0", 0, &status);
    else if (!strcmp(Hdr.gen.HdrVer,"Ver1.0.1"))
      fits_movnam_hdu(Fstokes, ASCII_TBL, "DUMPREF0", 0, &status);

    /* Get the current HDU number */
    fits_get_hdu_num(Fstokes, &NFirstTable);


    /* open new fits file for writing */
    /*    sprintf(RotFile,"%s.rot.fits", OutRoot);

	  if(fits_create_file(&Frot, RotFile, &fitsstatus)) {
	  printf("Error opening FITS file %s!!!\n",RotFile);fflush(stdout);
	  exit(4);
	  } */

    /* now start running through each dump */
    for(i_scan=0;i_scan<NDumps;i_scan++){

      /* move to next dump's data */
      if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
	fits_movabs_hdu(Fstokes, NFirstTable+i_scan, &hdutype, &status); 
      }
      else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
	fits_movabs_hdu(Fstokes,NFirstTable+(i_scan)*2+1,&hdutype,&status);
	fits_get_num_rows(Fstokes, &NPtsProf, &status);status=0; 
	fits_movrel_hdu(Fstokes, -1, NULL, &status);
      }


      /* find NPtsProf */
      ReadASPStokes(&Hdr, &StokesSubHdr, Fstokes, NPtsProf, 
		    Profiles, i_scan, Cmd->VerboseP);

      /* Check that StdBins equas NPtsProf from each input profile */
      if (NPtsProf != (long)StdBins) {
	printf("The number of bins in the data file %s (%ld) does ",
	       FitsFile, NPtsProf);
	printf("not equal that of the input standard profile (%d).  ",
	       StdBins);
	printf("Exiting...");
	exit(6);
      }

      for(i_chan=0;i_chan<Hdr.obs.NChan;i_chan++){

	/*********** READING DONE.  BEGIN PHASE DETERMINATION ************/

	  /* Bad scans are zeroed so if summ of the profile is zero, it's 
	     not to be used in summation */
	ProfSum = FSum(&Profiles[i_chan].rstds[0], NPtsProf);
	if(ProfSum != 0.0) { // i.e. good data
	  
	  
	  memcpy(profs,Profiles[i_chan].rstds,sizeof(float)*NBINMAX);
	  
	  /* Now fftfit to find shift required in second profile */
	  fftfit_(profs,&stdamps[1],&stdphas[1],
		  &StdBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);  
	  
	  /** Determine difference between read-in profile and standard 
	      profile **/
	  if (Cmd->VerboseP) 
	    printf("Shift, b: %f +/- %f, %f +/- %f\n\n",Shift,EShift,b,errb);
	  
	  /* Convert shift to phase (0,1) and radians (0,2pi) */
	  PhaseShift = (double)(-Shift/(double)StdBins);
	  AngleShift = (double)(-Shift/(double)StdBins*TWOPI);
	  
	  if(Cmd->VerboseP){
	    printf("Dump %d, Channel %d:\n",i_scan, i_chan);
	    printf("----------------------\n");
	    printf("Shift: %f +/- %f bins\n", Shift, EShift);
	    printf("          %lf +/- %lf rad\n",
		   AngleShift,EShift/StdBins*TWOPI);
	    printf("          %lf +/- %lf phase\n",PhaseShift, EShift/StdBins);
	    printf("Scale factor: %f +- %f \n\n",b,errb);
	  }
	  
	  /*********** PHASE DIFF COMPLETED ***********/
	  
	  /*** Calculate time to add to time stamp of scan:  ***/
	  
	  /* First make shift positive if it is not: */
	  if(Shift < 0.0) Shift +=StdBins;
	  
	 /* Change in rot. phase = (phase shift) - (initial phase of profile) */
	  PhaseChange = Shift/(double)StdBins-StokesSubHdr.DumpRefPhase[i_chan];
	  /* Make sure it's positive */
	  if (PhaseChange < 0.0) PhaseChange += 1.0;
	  
	  if(Cmd->VerboseP) 
	    printf("Initial phase: %lf\nPhase Change: %f +/- %f\n\n",
		   StokesSubHdr.DumpRefPhase[i_chan], 
		   PhaseChange, EShift/StdBins);
	  
	  /* Determine new fractional MJD: 
	     new time = old time + (pulse period)*(phase diff) */
	  MJD0 = Hdr.obs.IMJDStart;
	  FracMJD = (StokesSubHdr.DumpMiddleSecs + 
		     StokesSubHdr.DumpRefPeriod[i_chan]*PhaseChange)/86400.0;
	  /* Adjust integer MJD if FracMJD is < 0.0 or >= 1.0; also adjust 
	     FracMJD so it is between 0.0 and 1.0 */
	  if(FracMJD < 0.0){
	    MJD0 -= 1;
	    FracMJD += 1.0;
	  }
	  if(FracMJD >= 1.0){
	    MJD0 += 1;
	    FracMJD -= 1.0;
	  }
	  
	  TOAErr = StokesSubHdr.DumpRefPeriod[i_chan] *1.e6 *EShift/StdBins;
	  
	  //	MJDPaste(MJD0,FracMJD,TOA);
	  
	  /* Past together MJD into a string */
	  
	  if (Cmd->Tempo2P) { // can afford more decimal places :)
	    sprintf(frac,"%18.16lf",FracMJD);
	    sprintf(TOA,"%5d%17s",MJD0,&frac[1]);
	  }
	  else {
	    sprintf(frac,"%15.13lf",FracMJD);
	    sprintf(TOA,"%5d%14s",MJD0,&frac[1]);
	  }
	  NToa++; // count up TOAs
	  //	TotToa++;
	  
	  strncpy(Source,Hdr.target.PSRName,7);
	  strncpy(&Source[7],"\0",1);
	  
	  /* If writing in tempo2 format, set up flags */
	  if (Cmd->Tempo2P) {
	    strcpy(AllT2Flags,"\0");
	    if(Cmd->MJDFlagP) {
	      /* Make MJD flag argument have two decimal places */
	      strncpy(tempMJD, TOA, 8);
	      //	    strncpy(&tempMJD[8],"\0",1);
	      sprintf(AllT2Flags,"-mjd %8s",tempMJD);
	      
	    }
	    if(Cmd->BEFlagP) {
	      if (!strcmp(NObs,"1")) {
		sprintf(AllT2Flags, "%s -be gasp", AllT2Flags); 
	      } 
	      else if (!strcmp(NObs,"3")) {
		sprintf(AllT2Flags, "%s -be asp", AllT2Flags); 	      
	      }
	      else if (!strcmp(NObs,"f")) {
		sprintf(AllT2Flags, "%s -be bon", AllT2Flags); 	      
	      }
	      else {
		printf("Unrecognized observatory code %s. ", NObs);
		printf("Not implementing this flag\n");
	      }
	    }
	    /* Now stitch all the flags together, adding the extra 
	       user-definied flags if they exist */
	    if(Cmd->T2FlagsP){
	      sprintf(AllT2Flags, "%s %s", AllT2Flags, Cmd->T2Flags);
	    }
	  }
	  
	  freqwrite=mjdwrite=errwrite=0;
	  
	  /* If user did not restrict frequencies OR user did restrict 
	     frequencies AND frequency of current TOA is within restricted
	     range, then allow writing based on frequency */
	  if (!Cmd->FreqRangeP || ( Cmd->FreqRangeP && 
			   (Hdr.obs.ChanFreq[i_chan] >= Cmd->FreqRange[0] && 
			    Hdr.obs.ChanFreq[i_chan] <= Cmd->FreqRange[1]) ) )
	    freqwrite=1;
	  
	  if  (!Cmd->MJDRangeP || ( Cmd->MJDRangeP && 
				((double)MJD0+FracMJD >= Cmd->MJDRange[0] && 
				 (double)MJD0+FracMJD <= Cmd->MJDRange[1]) ) ) 
	    mjdwrite=1;

	  /* If user set a cutoff value for TOA error, then include the TOA 
	     only if the uncertainty is below the cutoff specified */
	  if  (!Cmd->ErrCutP || ( Cmd->ErrCutP && 
				  Cmd->ErrCut > TOAErr ) )
	    errwrite=1;

	  /*** write each new TOA to file in the correct tempo format ***/
	  
	  /* If it passes frequency, MJD, and toa uncertainty
	     restriction parameters then write to file */
	  if (freqwrite && mjdwrite && errwrite) {
	    NWrtToa++;
	    if (Cmd->Tempo2P) {
	      fprintf(Ftoa, "%s %8.3lf %22s %8.3f %3s %s\n",
		      Infile, Hdr.obs.ChanFreq[i_chan],
		      TOA, TOAErr, NObs, AllT2Flags);
	    }
	    else {
	      if (Cmd->NoIncrementP)
		fprintf(Ftoa, "%s%5d%8s%10.3f %19s%9.2f\n",
			NObs,1,Source,Hdr.obs.ChanFreq[i_chan],TOA,TOAErr);
	      else
		fprintf(Ftoa, "%s%5d%8s%10.3f %19s%9.2f\n",
			NObs,TotWrtToa+NWrtToa,Source,
			Hdr.obs.ChanFreq[i_chan],TOA,TOAErr);
	    }
	  }
	}
	else { // found an omitted scan
	  NOmit++;
	  if(Cmd->VerboseP || Cmd->CheckOmitP){
	    fprintf(Fcheck, "%6d     %.1lf\n",
		    i_scan,Hdr.obs.ChanFreq[i_chan]);
	  }
	}
      }
      
      if(Cmd->VerboseP) printf("\n");fflush(stdout);
    }
    
    /* Free memory used by Profiles structure since we will be re-malloc'ing
       it in the next loop-around */
    free(Profiles);
    
    
    //  for (i=0;i<Cmd->NFiles;i++)
    fits_close_file(Fstokes, &status);
    TotToa+=NToa;
    TotWrtToa+=NWrtToa;
    
    printf("%d TOAs found\n", NToa);
    if (Cmd->FreqRangeP || Cmd->MJDRangeP || Cmd->ErrCutP)
      printf("%d TOAs written based on command-line restrictions\n",NWrtToa);
    if (NOmit > 0) 
      printf("%d input scans were omitted from TOA calculation based on RFI filtering\n\n",NOmit);
    else
      printf("\n\n");
  }
  
  /******* END LOOP OVER FILES ********/
  if(Cmd->VerboseP || Cmd->CheckOmitP) fclose(Fcheck);
 
  /** close TOA file **/
  printf("+=========================================================================+\n\n");
  printf("Completed successfully.  Found %d TOAs.\n", TotToa);
  if (Cmd->FreqRangeP || Cmd->MJDRangeP || Cmd->ErrCutP)
    printf("Wrote %d to file based on command-line restrictions\n",TotWrtToa);
  printf("Output TOA file: %s \n\n",Toafile);fflush(stdout);
  
  fclose(Ftoa);
  
  
  exit(0);

}


/* void MJDPaste(int nmjd, double fmjd, char *toaout)
{
 
  char frac[16];
 
  sprintf(frac,"%15.13lf",fmjd);
  sprintf(toaout,"%5d%14s",nmjd,&frac[1]);
 
  } */
