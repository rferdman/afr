/* Scaled-down version of ASPFitsReader to simply read in the calibrated 
 * Stokes profiles, calculate Linear Polarization, and dump out an ascii file 
 * 
 * --R. Ferdman   August 6, 2004 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "ASPCommon.h"
#include "StokesCmdLine.h" 

int main(int argc, char *argv[])
{
  

  int StokesNDumps;
  int NFirstTable, NumHDU, hdutype;
  struct ASPHdr StokesHdr;
  struct SubHdr StokesSubHdr;

  struct StdProfs *StokesProfs;
  char StokesHead[256];
  int OutRootIndex=0, LastSlashIndex=0, FileNo;
  int DumpRange[2];
  char TempInFile[128],FitsFile[128], OutRoot[128], ProgName[32];
  char StokesFile[64], ParAngFile[64], OutFileEnd[16];
  fitsfile **Fstokes;
  FILE *fpStokes;
  FILE *fpParAng[NCHMAX];

  int i,j,nscan,status=0;
  int fitsstatus=0;
  long NPtsProf=0;
  double ParAng, MJD, RA, Dec;
  char NObs[3];

  double SBase,Srms,VBase,Vrms,LinAvg,Linrms,UBase,QBase,Urms,Qrms;
  double Duty,SPeak,p0218,LinPeak;
  float s0218[NBINMAX];
  double FinalMask[NBINMAX];
  int spk, RefBin;
  double x,ptype;
  Cmdline         *StokesCmd;

  /* Get command line variables */
  StokesCmd = parseCmdline(argc, argv);

  /* Normally use this somewhere, and not showOptionValues */
  StokesCmd->tool = StokesCmd->tool;
                                                                             
  /* Store program name */
  strcpy(ProgName, argv[0]);

  strcpy(FitsFile, StokesCmd->Infile);
  
  LastSlashIndex = -1;
  for(i=0;i<strlen(FitsFile);i++){
    if(!strncmp(&FitsFile[i],"/",1))
      LastSlashIndex = i;
  }
  /* Get output file root name from input file name */
  for (i=0;i<strlen(FitsFile);i++){
    if(!(strcmp(&FitsFile[i],"stokes.fits"))  ){
      OutRootIndex = i-1;
      sprintf(OutFileEnd,"stokes.asc");
      break;
    }
    if(!(strcmp(&FitsFile[i],"rot.fits")) ) {
      OutRootIndex = i-1;
       sprintf(OutFileEnd,"rot.asc");

      break;
    }
  }

  if(OutRootIndex == 0){
    printf("WARNING: input file name does not follow .stokes.fits or .rot.fits convention.\n");
    sprintf(OutRoot,"ASPOut");
  }
  else {
    strcpy(OutRoot,"\0");
    strncpy(OutRoot,&FitsFile[LastSlashIndex+1],OutRootIndex-LastSlashIndex-1);
  }

  /* Malloc fitsfile descriptor to number of data files */
  Fstokes  = (fitsfile **)malloc(StokesCmd->NFiles*sizeof(fitsfile));

  /* Open fits data file(s) */

  for(i=0;i<StokesCmd->NFiles;i++){
    status=0;
    if(StokesCmd->NFiles > 1){
      strncpy(TempInFile,&FitsFile[0],OutRootIndex);
      strcpy(FitsFile,"\0");
      sprintf(FitsFile,"%s.%d.stokes.fits",TempInFile,i);
    }

    if(fits_open_file(&Fstokes[i], FitsFile, READONLY, &status)){
      printf("Error opening FITS file %s !!!\n", FitsFile);
      exit(1);
    }
  }



  /* Read in values for header variables */
  if(ReadASPHdr(&StokesHdr, Fstokes[0]) < 0){
    printf("%s> Unable to read Header from file %s.  Exiting...\n",
	   ProgName,FitsFile);
    exit(1);
  }

  printf("\n==========================\n");
  printf("ASP FITS Header %s\n",StokesHdr.gen.HdrVer);
  printf("==========================\n\n");fflush(stdout);
  
  printf("Input file:  %s\n\n",FitsFile);fflush(stdout);

  printf("PSR %s:\n",StokesHdr.target.PSRName);
  printf("--------------\n\n");
  printf("Centre Frequency: %6.1lf MHz\n\n",StokesHdr.obs.FSkyCent);
  fflush(stdout);
  if(StokesCmd->NoBaseP)
    printf("Baseline subtraction turned OFF.\n\n");fflush(stdout);

  
  StokesNDumps = 0;
  for(i=0;i<StokesCmd->NFiles;i++){
    /* Get number of HDUs in fits file */
    fits_get_num_hdus(Fstokes[i], &NumHDU, &status);

    if(!strcmp(StokesHdr.gen.HdrVer,"Ver1.0")){
      StokesNDumps += NumHDU-3;  /* the "3" is temporary, depending on how 
				   many non-data tables we will be using */
    }
    else if(!strcmp(StokesHdr.gen.HdrVer,"Ver1.0.1")){
      StokesNDumps += (NumHDU-3)/2;
    }
    else{
      printf("Do not recognize FITS file version number in header.\n");
      printf("This header %s. Exiting...\n",StokesHdr.gen.HdrVer);
      fflush(stdout);
      exit(3);
    }
  }
 
  printf("Number of channels:  %d\n",StokesHdr.obs.NChan);
  printf("Number of dumps:     %d\n\n",StokesNDumps);

  /* Set dump ranges to output as ascii files */
  if (StokesCmd->DumpRangeP){
    DumpRange[0] = StokesCmd->DumpRange[0];
    DumpRange[1] = StokesCmd->DumpRange[1];
    printf("Will only output ascii files for dumps %d to %d.\n\n",
	   DumpRange[0],DumpRange[1]);
  }
  else{ // output all dumps
    DumpRange[0] = 0;
    DumpRange[1] = StokesNDumps-1;
  }



  StokesProfs = (struct StdProfs *)malloc
    (StokesHdr.obs.NChan*sizeof(struct StdProfs));

  /* Convert Ra and Dec, mainly for use with parallactic angle calculation if needed */
  if(StokesHdr.target.RA > 0. && fabs(StokesHdr.target.Dec) > 0. ){
    RA = (StokesHdr.target.RA/24.)*TWOPI;
    Dec = StokesHdr.target.Dec*TWOPI/360.;
  }
  else{
    printf("\nRA and Dec are not in file. Parallactic angle calculations will be\n");
    printf("incorrect...\n");fflush(stdout);
  }

  /* Get observatory code */
  sscanf(StokesHdr.obs.ObsvtyCode,"%s",NObs);

  /* Move to the first data table HDU in the fits file */
  if(!strcmp(StokesHdr.gen.HdrVer,"Ver1.0"))
    fits_movnam_hdu(Fstokes[0], BINARY_TBL, "STOKES0", 0, &status);
  else if (!strcmp(StokesHdr.gen.HdrVer,"Ver1.0.1"))
    fits_movnam_hdu(Fstokes[0], ASCII_TBL, "DUMPREF0", 0, &status);

  /* Get the current HDU number */
  fits_get_hdu_num(Fstokes[0], &NFirstTable);

  FileNo = 0;
  for(nscan=0;nscan<StokesNDumps;nscan++){

    if (nscan>=DumpRange[0] && nscan<=DumpRange[1]) {
      /* move to next dump's data */
      if(!strcmp(StokesHdr.gen.HdrVer,"Ver1.0")){
	fits_movabs_hdu(Fstokes[FileNo], NFirstTable+nscan, &hdutype, &status); 
      }
      else if(!strcmp(StokesHdr.gen.HdrVer,"Ver1.0.1")){
	/* if we've reached the end of the FITS file then increase FileNo */
	if( (nscan%MAXDUMPS == 0) && (nscan > 0) ) {
	  FileNo++;
	  printf("Switching to File %d\n",FileNo);fflush(stdout);
	}
	fits_movabs_hdu(Fstokes[FileNo],NFirstTable+(nscan%MAXDUMPS)*2+1,
			&hdutype,&status);
	fits_get_num_rows(Fstokes[FileNo], &NPtsProf, &status);status=0; 
	fits_movrel_hdu(Fstokes[FileNo], -1, NULL, &status);
      }


      /* find NPtsProf */
      ReadASPStokes(&StokesHdr, &StokesSubHdr, Fstokes[FileNo], NPtsProf, 
		    StokesProfs, nscan, StokesCmd->VerboseP);

      /* pspmtoa header format */
      for(i=0;i<StokesHdr.obs.NChan;i++){

	if(StokesCmd->HeaderP){
	  sprintf(StokesHead,"%4.0f. %4d   %7s            1.000",
		  StokesHdr.obs.FSkyCent,StokesHdr.redn.RNBinTimeDump, 
		  StokesHdr.target.PSRName);	
	}
	else {
	  sprintf(StokesHead,"# %.1f %.7f %.10f %ld %.3f %.3f %d %s %d %9s %.10f",
		  (double)StokesHdr.obs.IMJDStart, StokesSubHdr.DumpMiddleSecs, 
		  StokesSubHdr.DumpRefPeriod[i],(long)1,
		  StokesHdr.obs.ChanFreq[i], 
		  StokesHdr.obs.DM, StokesHdr.redn.RNBinTimeDump,
		  StokesHdr.obs.ObsvtyCode, 1, StokesHdr.target.PSRName, 
		  StokesSubHdr.DumpRefPhase[i]);             
	}

	sprintf(StokesFile,"%s.%4.4d.%4.4d.%4.4d.%s",
		OutRoot,(int)StokesHdr.obs.ChanFreq[i],
		(int)(10000*(StokesHdr.obs.ChanFreq[i]
			     -(int)StokesHdr.obs.ChanFreq[i])),
		nscan,OutFileEnd);
	if ((fpStokes = fopen(StokesFile,"w")) == 0)
	  { printf("Cannot open %s. Exiting...\n",StokesFile); exit(1); }
	fprintf(fpStokes,"%s\n",StokesHead);
      
	Duty = DutyLookup(StokesHdr.target.PSRName);
	BMask(StokesProfs[i].rstds,&StokesHdr.redn.RNBinTimeDump,&Duty,FinalMask);
	Baseline(StokesProfs[i].rstds,FinalMask,&StokesHdr.redn.RNBinTimeDump,
		 &SBase,&Srms);
	Baseline(StokesProfs[i].rstdv,FinalMask,&StokesHdr.redn.RNBinTimeDump,
		 &VBase,&Vrms);
	Baseline(StokesProfs[i].rstdq,FinalMask,&StokesHdr.redn.RNBinTimeDump,
		 &QBase,&Qrms);
	Baseline(StokesProfs[i].rstdu,FinalMask,&StokesHdr.redn.RNBinTimeDump,
		 &UBase,&Urms);
      
	/* special provision */
	if(strncmp(StokesHdr.target.PSRName,"0218",4) == 0) {  
	  printf("Changing 0218 baseline.\n");
	  for(j=0;j<StokesHdr.redn.RNBinTimeDump;j++) 
	    s0218[j] = -StokesProfs[i].rstds[j];
	  p0218 = FindPeak(s0218,&StokesHdr.redn.RNBinTimeDump,&spk);
	  SBase = -p0218;
	}
      
	if(!StokesCmd->NoBaseP){
	  for(j=0;j<StokesHdr.redn.RNBinTimeDump;j++) {
	    StokesProfs[i].rstds[j] -= SBase;
	    StokesProfs[i].rstdq[j] -= QBase;
	    StokesProfs[i].rstdu[j] -= UBase;
	    StokesProfs[i].rstdv[j] -= VBase;  
	  }
	}      

	SPeak =  FindPeak(StokesProfs[i].rstds,&StokesHdr.redn.RNBinTimeDump,
			  &spk);
	StokesProfs[i].SNR = SPeak*Srms;

	/* Calculate linear polarization, position angle, and PA error */

	if(i==0){
	  printf("Writing ascii files %s.*.%4.4d.%s\n\n",OutRoot,nscan,
		 OutFileEnd);
	  fflush(stdout);
	}

	for(j=0;j<StokesHdr.redn.RNBinTimeDump;j++) {
	  /* Q^2 + U^2 */
	  StokesProfs[i].stdlin[j] = sqrt(StokesProfs[i].rstdq[j]*
					  StokesProfs[i].rstdq[j]+ 
					  StokesProfs[i].rstdu[j]*
					  StokesProfs[i].rstdu[j]);  
	  /* arctan(U/Q) -- factor of two taken into account for file writing */
	  StokesProfs[i].stdphi[j] = -atan2(StokesProfs[i].rstdu[j],
					    StokesProfs[i].rstdq[j]);  

	  StokesProfs[i].stdphierr[j] = sqrt((StokesProfs[i].rstdq[j]*
					      StokesProfs[i].rstdq[j]
					      /(Urms*Urms) +
					      StokesProfs[i].rstdu[j]*
					      StokesProfs[i].rstdu[j]
					      /(Qrms*Qrms))/
					     (StokesProfs[i].stdlin[j]*
					      StokesProfs[i].stdlin[j]*
					      StokesProfs[i].stdlin[j]*
					      StokesProfs[i].stdlin[j]));  
	}
	Baseline(StokesProfs[i].stdlin,FinalMask,&StokesHdr.redn.RNBinTimeDump,
		 &LinAvg,&Linrms);
	for(j=0;j<StokesHdr.redn.RNBinTimeDump;j++) {
	  StokesProfs[i].stdlin[j] -= LinAvg;
	}
      

	if(StokesCmd->VerboseP){
	  printf("Q, U, lin, total rms: %f %f %f %f\n",1./Qrms, 1./Urms, 
		 1./Linrms,1./Srms);
	  printf("SNR: %f\n",StokesProfs[i].SNR);
	}

	for(j=0;j<StokesHdr.redn.RNBinTimeDump;j++) {
	  /* see how strong the linear polarization is */
	  x = StokesProfs[i].stdlin[j]*Srms; 
	  ptype = 43.1;
	  if (x > 1.) ptype=43.2;
	  if (x > 2.) ptype=43.3;
	  if (x > 3.) ptype=43.4;
	  if (x > 4.) ptype=43.5;
	  if (x > 5.) ptype=43.6;
	  fprintf(fpStokes,
		  "%5d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%6.1f\n",j,
		  StokesProfs[i].rstds[j],StokesProfs[i].rstdq[j],
		  StokesProfs[i].rstdu[j],StokesProfs[i].rstdv[j],
		  /* phi in degrees */
		  StokesProfs[i].stdlin[j],
		  StokesProfs[i].stdphi[j]*180.0/TWOPI, 
		  StokesProfs[i].stdphierr[j]*180.0/TWOPI,ptype);
	}
	fclose(fpStokes);

	if(StokesCmd->ParAngFileP){
	  /* open file for writing -- one for each channel */
	  if(nscan==0) {
	    sprintf(ParAngFile,"%s.%4.4d.%4.4d.parang.asc",
		    OutRoot,(int)StokesHdr.obs.ChanFreq[i],
		    (int)(10000*(StokesHdr.obs.ChanFreq[i]-
				 (int)StokesHdr.obs.ChanFreq[i])));
	    /* take peak S of this first scan as reference bin */
	    LinPeak = FindPeak(StokesProfs[i].stdlin,
			       &StokesHdr.redn.RNBinTimeDump,&RefBin);
	    printf("\nCreating parallactic angle vs. Stokes parameters files\n");
	    printf("Reference bin (of peak S): %d\n",spk);
	    if ((fpParAng[i] = fopen(ParAngFile,"w")) == 0)
	      { printf("Cannot open %s. Exiting...\n",ParAngFile); exit(1); }
	    /* Make header for parallactic angle files */
	    fprintf(fpParAng[i],"%s  %7.2lf  %4.4d  %8.2lf\n",
		    StokesHdr.target.PSRName,StokesHdr.obs.ChanFreq[i],
		    RefBin,LinPeak); 
	  }
	  /* Now calculate parallactic angle */
	  MJD = (double)StokesHdr.obs.IMJDStart + 
	    (StokesSubHdr.DumpMiddleSecs)/86400.;
	  ParAng = GetChi(StokesHdr.target.PSRName,MJD,NObs,RA,Dec);
	  fprintf(fpParAng[i],"%14.5lf%15.7f%15.7f%15.7f%15.7f%15.7f\n", 
		  MJD, ParAng*180.0/TWOPI,
		  StokesProfs[i].rstds[RefBin],StokesProfs[i].rstdq[RefBin],
		  StokesProfs[i].rstdu[RefBin],StokesProfs[i].rstdv[RefBin]);
	}


      }

      


    if(StokesCmd->VerboseP) printf("\n");fflush(stdout);

  }// DumpRange if condition
  }


  for (i=0;i<StokesCmd->NFiles;i++)
    fits_close_file(Fstokes[i], &fitsstatus);

    /* Close parallactic angle file */
  if(StokesCmd->ParAngFileP)
    for(j=0;j<StokesHdr.obs.NChan;j++)
      fclose(fpParAng[j]);

  printf("Completed successfully.\n\n");fflush(stdout);

  exit(0);

}
