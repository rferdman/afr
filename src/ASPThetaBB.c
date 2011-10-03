/********************************************************************/
/*                                                                  */
/*  ThetaBB:                                                        */
/*                                                                  */
/*  - Will calculate ThetaYY to correct phase offset between        */
/*    polarizations                                                 */
/*                                                                  */
/* Robert Ferdman, University of British Columbia, July 19, 2005    */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "ASPCommon.h"
#include "ThetaCmdLine.h" 
#include "ASPCal.h"


int main(int argc, char **argv)
{
  
  int            i, j, status=0, NumHDU;

  double         *ASquared[NCHMAX], *BSquared[NCHMAX];
  double         *ReAconjB[NCHMAX], *ImAconjB[NCHMAX];

  struct ASPHdr  Hdr;
  struct SubHdr  *SubHdr;
  struct RunVars RunMode;
  struct CalVars CalMode;


//  int            PhaseBin[2];
  int            OnBin[2], OffBin[2];

  double         OnAvgRe, OnAvgIm;
  double         OffAvgRe, OffAvgIm;
  double         CalHeightRe, CalHeightIm;

  double         ThetaBB[NCHMAX];
  char           ProgName[256];

  fitsfile       *Fin;
  char           Outfile[32];
  FILE           *Fout;
  Cmdline         *Cmd;

  float          ThetaDeg[NCHMAX], Chan[NCHMAX];

  /* Variables associated with 'pacv -b' output file */
  int            i_b, ThisChan;
  double         Bfreq, C1, C2, C3, C1err, C2err, C3err, Bangle[NCHMAX];
  char           BLine[256];
  FILE           *Fb;

  int Diagnose=1;
  FILE *fptest;

  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  
  /*   showOptionValues();  */

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

                                       
  /* Make sure that there is one and only one argument */
  /*  if (argc != 2){
    printf("\nUsage: ThetaBB <input file>\n\n");fflush(stdout); 
    exit(2);
    } */
                            
  /* Store program name */
  strcpy(ProgName,argv[0]); 
  /* Store input file name */ 
  //  strcpy(RunMode.Infile,argv[1]);

  /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  strcpy(RunMode.Infile,Cmd->Infile);


  /* Fill up some necessary RunMode variables */
  RunMode.AddDumps = 0;
  RunMode.Verbose = 1;
  RunMode.OldFits = 0;

  /* Initialize header variables */
  InitPars(&Hdr);

  /* Grab pulsar cal file name */
  //  strcpy(RunMode.Infile,CalCmdCalfile);

  /* Open pulsar cal fits file */
  status=0;

  if(fits_open_file(&Fin, RunMode.Infile, READONLY, &status)){
    printf("Error opening FITS cal file %s !!!\n", RunMode.Infile);
    exit(1);
  }

  /* Read in values for header variables */
  if(ReadASPHdr(&Hdr, Fin) < 0){
    printf("Unable to read Header from CAL file %s.  Exiting...\n",
	   RunMode.Infile);
    exit(2);
  };

  /* Figure out sideband */
  if (Hdr.obs.ChanFreq[0] < Hdr.obs.ChanFreq[Hdr.obs.NChan-1]){
    RunMode.Sideband = 1;   /* upper sideband */
    printf("\n\n Data is UPPER sideband.\n\n");fflush(stdout);
  }
  else{
    RunMode.Sideband = -1;  /* lower sideband */
    printf("\n\n Data is LOWER sideband.\n\n");fflush(stdout);
  }

  /* Get number of HDUs in fits file */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  RunMode.NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				    many non-data tables we will be using */
  if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
    RunMode.NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				   many non-data tables we will be using */
  }
  else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
    RunMode.NDumps = (NumHDU-3)/2;
    printf("Number of dumps in file %s:  %d\n",RunMode.Infile,
	   (NumHDU-3)/2);
  }
  else{
    printf("Do not recognize FITS file version number in header.\n");
    printf("This header %s. Exiting...\n",Hdr.gen.HdrVer);fflush(stdout);
    exit(3);
  }
  
  /* Make sure that this *is* a cal scan! */
  if(strcmp(Hdr.gen.ObsMode,"CAL") != 0){
    printf("File %s is NOT a cal scan according to header!  Exiting...\n",
           RunMode.Infile);fflush(stdout);
    return -2;
  }

  RunMode.NBins = Hdr.redn.RNBinTimeDump;

  /* Read in cal data */
  if(GetCalData(&Hdr, SubHdr, &RunMode, Fin, 
		&ASquared[0], &BSquared[0], &ReAconjB[0], &ImAconjB[0]) < 0) {
    printf("Error: Could not read cal data from file %s.  Exiting...\n",
	   RunMode.Infile);fflush(stdout);
    exit(4);
  }

  /* Open output file which will be read into ASPFitsReader */
  sprintf(Outfile,"thetaBB.out");
  if((Fout = fopen(Outfile,"w")) == NULL){
    printf("Could not open file %s.  Exiting...\n",Outfile);fflush(stdout);
    exit(6);
  }

  
  /* Find ranges of bins for on and off for each channel -- default is all */
  if(GetPhases(&Hdr, &RunMode, &CalMode, 
		&ASquared[0], &BSquared[0], OnBin, OffBin) <0){
    printf("Problem getting on/off phases for cal file %s. Exiting...\n",
	   RunMode.Infile);fflush(stdout);
    exit(5);
  }
  

  /* Read in output file from pacv -b */
  if(Cmd->BfileP){
    printf("Reading input file from 'pacv -b' output, named %s\n",
	   Cmd->Bfile);
    if((Fb = fopen(Cmd->Bfile, "r")) == NULL){
      printf("Cannot open file %s for reading in 'pacv -b' output.  Exiting...\n",
	     Cmd->Bfile);fflush(stdout);
      return -2;
    }
    printf("File %s open.\n", Cmd->Bfile);
    while (fgets(BLine, 256, Fb) != NULL){
      sscanf(BLine, "%d%lf%lf%lf%lf%lf%lf%lf", 
	     &i_b, &Bfreq, &C1, &C1err, &C2, &C2err, &C3, &C3err);

      /* Now do frequency match to determine corresponding channel index */
      if((ThisChan = Freq2Chan(Bfreq, Hdr.obs.ChanFreq, Hdr.obs.NChan)) < 0) {
	/* If mapping unsuccessful, quit out with error */
	fprintf(stderr, "Error in 'pacv -b' output. Frequency %f does not ", 
		Bfreq);
	fprintf(stderr, "map to any channel for this cal data file. \n");
	return -3;
      }


      /* Assuming here feeds are linear?  Or is this general? If so,
       are we then assuming linearity for the ThetaBB calculation, which uses 
       Re and Im cross products without checking polarization basis? */
      /* C2 = U; C3 = V */
      Bangle[ThisChan] = atan2(C2, C3);
      //printf("Bangle[%d] = %lf\n", ThisChan, Bangle[ThisChan]);
	
    }

  }
  fclose(Fb);
  printf("File %s closed.\n", Cmd->Bfile);


  /*** BEGIN FIND THETA_BB ***/

  printf("\n\nResults from ThetaBB calculation:\n");
  printf("-------------------------------------\n\n");
  printf(" Freq(MHz)  ThetaBB(deg)  Error on ThetaBB(deg)\n");
  printf(" ---------  ------------  ---------------------\n");


  if(Diagnose){
    if((fptest = fopen("ThetaBB_CalHeight.dat","w")) == NULL){
      printf("Could not open file ThetaBB_CalHeight.dat. Exiting...\n");
      fflush(stdout);
      exit(1);
    }
    fprintf(fptest,"Freq      CalHeightRe    OnAvgRe       OffAvgRe    CalHeightIm    OnAvgIm        OffAvgIm\n\n");
  }
  

  
  for(j=0;j<Hdr.obs.NChan;j++){

    /* if lower sideband, do Im(L*R) --> -Im(L*R)*/
    /* if SWAP --> Im(L*R) --> -Im(L*R) */
    for(i=0;i<RunMode.NBins;i++) {
      ImAconjB[j][i] *= RunMode.Sideband;
      if (Cmd->SwapP)
	ImAconjB[j][i] *= -1.0;
    }

    CalHeightRe = GetCalHeight(ReAconjB[j], 
			       RunMode.NBins,OnBin,OffBin,&OnAvgRe,&OffAvgRe);

    CalHeightIm = GetCalHeight(ImAconjB[j], 
			       RunMode.NBins,OnBin,OffBin,&OnAvgIm,&OffAvgIm);

    if(Diagnose){
      fprintf(fptest,"%6.1lf    %11.3f  %11.3f  %11.3f    %11.3f  %11.3f   %11.3f\n",
	      Hdr.obs.ChanFreq[j],CalHeightRe,OnAvgRe,OffAvgRe,
	      CalHeightIm,OnAvgIm,OffAvgIm);
    }
      

    ThetaBB[j] = atan2(CalHeightIm, CalHeightRe); // * 360./TWOPI;

    /* Remove extra angle if using output file from "pacv -b" */
    if(Cmd->BfileP)
      ThetaBB[j] -= Bangle[j];

    /* do error calculations -- need RMS as well -- see Shauna's thesis */

    ThetaDeg[j] = (float)(ThetaBB[j]*360./TWOPI);
    Chan[j] = (float)Hdr.obs.ChanFreq[j];
  
    printf("%10.1lf  %12.2f\n",Hdr.obs.ChanFreq[j],ThetaDeg[j]);
    fflush(stdout);
    fprintf(Fout,"%8.1lf %8.5lf\n",Hdr.obs.ChanFreq[j],ThetaBB[j]);
    
  }

  if(Diagnose) fclose(fptest);
  
  fits_close_file(Fin, &status);
  fclose(Fout);

    
exit(0);

}


