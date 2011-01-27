/* Takes output stokes profile fits file and rotate it by an angle provided 
 * by user 
 * 
 * --R. Ferdman,  13 May 2008 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "RotateCmdLine.h" 
#include "ASPCommon.h"

int main(int argc, char *argv[])
{
  

  int NDumps;
  int NFirstTable, NumHDU, hdutype;
  struct ASPHdr Hdr;
  struct SubHdr StokesSubHdr;

  struct StdProfs *Profile;
  int OutRootIndex=0, LastSlashIndex=0, TotalRootIndex=0;
  char FitsFile[128], OutRoot[128], ProgName[32];
  char RotFile[64];
  fitsfile *Fstokes;
  fitsfile *Frot;

  int i,nscan,status=0;
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
  sprintf(RotFile,"%s.rot.fits", OutRoot);

  if(fits_create_file(&Frot, RotFile, &fitsstatus)) {
    printf("Error opening FITS file %s!!!\n",RotFile);fflush(stdout);
    exit(4);
  }

  printf("Rotating profile(s) by %6.3lf radians or %6.2lf degrees...\n\n",
	 Cmd->ByAngle,Cmd->ByAngle*360./TWOPI);

  printf("Writing rotated file %s\n\n",RotFile);fflush(stdout);

  /* Write header info into file */
  if(WrtASPHdr(&Hdr, Frot) < 0) {
    printf("Error writing ASP header structure!\n");fflush(stdout);
    exit(5);
  }
  
   /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  strcpy(RunMode.Infile,Cmd->Infile); 

 /* Set up RunMode structure to be compatible with WrtASPStokes */
  RunMode.AddChans = 0;
  RunMode.AddDumps = 0;
  RunMode.NScanOmit = 0;
  IZero(RunMode.DumpOmit,MAXOMIT);
  IZero(RunMode.ChanOmit,MAXOMIT);
  IZero(RunMode.ZapChan,Hdr.obs.NChan);
  RunMode.NBins = RunMode.NBinsOut = Hdr.redn.RNBinTimeDump;
  
  
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

      /*********** READING DONE.  BEGIN ROTATION ************/
      
      RotateProf(&RunMode, &Profile[i], Cmd->ByAngle);
      
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
	
    }

    if(WrtASPStokes(Hdr, StokesSubHdr, Frot, nscan, 
		    Profile, &RunMode) < 0) {
      printf("Cannot write data tables to fits file. Exiting...\n");fflush(stdout);
      exit(6);
    }
    
   if(Cmd->VerboseP) printf("\n");fflush(stdout);
  }


  //  for (i=0;i<Cmd->NFiles;i++)
  fits_close_file(Fstokes, &fitsstatus);
  fits_close_file(Frot, &fitsstatus);

  printf("Completed successfully.\n\n");fflush(stdout);

  exit(0);

}
