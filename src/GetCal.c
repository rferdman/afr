#include <stdio.h>
#include <math.h>
/* #include "ASPHeader.h" */
#include "fitsio.h"
#include "ASPCommon.h"

void GetCalFac(struct ASPHdr *hdr, struct RunVars *RunMode, double *Profile, 
	       int NPtsProf, double CalValue, double *CalFactor, double *Tsys, 
	       double *Gain)
{
  int j;
  double Diff;
  /*   double ByPhase; */
  double MaxDiff1, MaxDiff2;
  int UpBin, DownBin;
  int PulseInside;
  int NLow, NHigh;
  double LowAvg, HighAvg;
  double CountsPerCal;
  double CalHeight;
  
  

  /* First rotate cal profile to centre on 180 degrees */
  /*   cprofc(Profile, NPtsProf, Amp, Phase);  */
  /*   ByPhase = Phase[1] - TWOPI/2.;  */
  /*   RotateProf(RunMode, Profile, ByPhase);      */


  /* first find largest differences between count values -- 
     these are the bins to mark */
  MaxDiff1 = 0.;
  MaxDiff2 = 0.;
  UpBin = 0;
  DownBin = 0;
      
  for (j=0;j<NPtsProf;j++) {
    if (j==0)
      Diff = Profile[0] - Profile[NPtsProf-1];
    else
      Diff = Profile[j] - Profile[j-1];
    /*     printf("%d: Diff = %f, fabs(Diff) = %f, MaxDiff1 = %f, 
	   MaxDiff2 = %f\n",j, Diff, fabsf(Diff), MaxDiff1, MaxDiff2); */
    if(fabsf(Diff) > fabs(MaxDiff1)){ 
      MaxDiff2 = MaxDiff1;
      MaxDiff1 = Diff;
      
      if (MaxDiff1 > 0.){ 
	UpBin = j;
      }
      else{
	DownBin = j;
      }
    }

    else if (fabsf(Diff) > fabs(MaxDiff2)){
      MaxDiff2 = Diff;
      if (MaxDiff2 > 0.){ 
	UpBin = j;
      }
      else{
	DownBin = j;
      }
    }

  }
  /*   printf("MaxDiff1 = %f, MaxDiff2 = %f\n\n",
       MaxDiff1, MaxDiff2);fflush(stdout); */

  /* stdin entry */
  /* printf("Enter UpBin for this freq: \n");fflush(stdout);
  fscanf(stdin,"%d",&UpBin);
  printf("Enter DownBin for this freq: \n");fflush(stdout);
  fscanf(stdin,"%d",&DownBin); */


  printf("GETCALFAC: UpBin = %d, DownBin = %d\n",UpBin, DownBin);fflush(stdout);
  
  NLow  = 0;
  NHigh = 0;
  LowAvg  = 0.;
  HighAvg = 0.;

  PulseInside = (DownBin > UpBin);  
  /* this refers to if cal pulse is entirely inside the profile 
   *
   * ie.    PulseInside == 1       PulseInside == 0
   *              ____              ____      ____
   *             |    |                 |    |
   *             |    |                 |    |
   *         ____|    |____             |____|
   *
   *
   *  This way, we test first if index is within the crest (if 
   *  PulseInside==1) or trough (if PulseInside==0).
   *
   */

  for (j=0;j<NPtsProf;j++){
    if(PulseInside){

      if(j>UpBin+1 && j<DownBin-1){  /* Top of pulse */
	NHigh++;
	HighAvg += Profile[j];
      }
      else if(j<UpBin-1 || j>DownBin+1){
	NLow++;
	LowAvg  += Profile[j];	
      }

    }
    else{

      if(j>DownBin+1 && j<UpBin-1){  /* Bottom of pulse */
	NLow++;
	LowAvg  += Profile[j];	
      } 
      else if(j<DownBin-1 || j>UpBin+1){
	NHigh++;
	HighAvg += Profile[j];
      }

    }
    
  }
  LowAvg  /= NLow;
  HighAvg /= NHigh;
  
  CalHeight = HighAvg - LowAvg; 

  /* Figure out Tsys and counts per cal given user-provided tcal
   * this is for GB only for now.

   * FUTURE: implement Don's engineering tcal extractor for flexibility -- 
   *         then won't need to specify tcals on command line.
   * ALSO:   do this for Arecibo!
   * ALSO:   assuming gain now, but maybe implement calculation depending on 
   *         elevation angle  */

  if (!strcmp(hdr->obs.ObsvtyCode,"1"))    /* GBT */
    *Gain = 2.;   /*  K/Jy; assuming 2 for now for GBT */
  /* For Arecibo -- need more sophisticated way of doing this...
   * these are averages right now */
  else if(!strcmp(hdr->obs.ObsvtyCode,"3")){
    if(!strcmp(hdr->gen.FEName, "327") )        *Gain = 11.;
    else if(!strcmp(hdr->gen.FEName, "430") )   *Gain = 11.;
    /* carriage house?  gain = 20.  */
    else if(!strcmp(hdr->gen.FEName, "430CH") ) *Gain = 20.;
    else if(!strcmp(hdr->gen.FEName, "610") )   *Gain = 11.;
    else if(!strcmp(hdr->gen.FEName, "ALFA") )  *Gain = 10.;
    else if(!strcmp(hdr->gen.FEName, "LWIDE") ) *Gain = 10.;
    else if(!strcmp(hdr->gen.FEName, "SLOW") )  *Gain = 8.;
    else if(!strcmp(hdr->gen.FEName, "SHIGH") ) *Gain = 8.5;
    else {
      printf("Receiver %s not in program, so don't know gain.  Exiting...\n",
	     hdr->gen.FEName);fflush(stdout);
      exit(1); 
    }
    

  }


  if(RunMode->CalMethod==0){
    *CalFactor = CalValue / CalHeight;
  }
/*   else if(RunMode->CalMethod==1) { */
  else {
    /*     CountsPerCal = CalValue; */
    *Tsys = (LowAvg/CalHeight)*CalValue; 
    *CalFactor = CalValue/((*Gain)*CalHeight); 
  }

  if(RunMOde->Verbose)
    printf("GETCALFAC:  CalHeight = %lf\n",CalHeight);fflush(stdout);

  /*   printf("GETCALFAC:  HighAvg = %f, LowAvg = %f\n",HighAvg,LowAvg); */


}





void GetCal(struct ASPHdr *hdr, struct RunVars *RunMode, double **CalFactor)
{

  int     i, j, status=0, colnum, anynull;
  int     NColumns;
  /*  long    NPtsProf=0; */
  /*   float  *LSquared, *RSquared, *ReLconjR, *ImLconjR; */
  struct ASPHdr CalHdr;
  long  NPtsProf=0;
  char Calfile[128];
  double **LSquared, **RSquared, **ReLconjR, **ImLconjR;
  int **SampleCount;
  fitsfile *Fcal;
  double Tsys, Gain;




  if(!RunMode->Cal){  /* if we are not calibrating */
    for (i=0;i<hdr->obs.NChan;i++)
      CalFactor[i][0] =CalFactor[i][1] =CalFactor[i][2] =CalFactor[i][3] = 1.0;
    return;
  }

  printf("GETCAL:  CalMethod = %d\n",RunMode->CalMethod);fflush(stdout);

  strcpy(Calfile,  RunMode->Calfile);

  /* Open cal fits file */
  if(fits_open_file(&Fcal, Calfile, READONLY, &status)){
    printf("Error opening FITS cal file %s !!!\n", RunMode->Calfile);
    exit(0);
  }

  if(ReadASPHdr(&CalHdr, Fcal) < 0){
    printf("Unable to read Header from CAL file %s.  Exiting...\n",Calfile);
    exit(1);
  };

  if(strcmp(CalHdr.gen.ObsMode,"CAL") != 0){
    printf("File %s is NOT a cal scan according to header!  Exiting...\n",
	   Calfile);fflush(stdout);
    exit(1);
  }


  fits_movnam_hdu(Fcal, BINARY_TBL, "ASPOUT0", 0, &status);   
  fits_get_num_rows(Fcal, &NPtsProf, &status);status=0; /* find NPtsProf */
/*   RunMode->NBins = (int)NPtsProf; */
  fits_get_num_cols(Fcal, &NColumns, &status);status=0;

  printf("CALFILE:  NChan = %d, NPtsProf = %ld, NColumns = %d\n",
	 hdr->obs.NChan,NPtsProf,NColumns);
  printf("          CalValue[0] = %f, CalValue[1] = %f\n",
	 RunMode->CalValue[0],RunMode->CalValue[1]);fflush(stdout); 

  /* malloc array for each channel based on profile size */

  LSquared = (double **)malloc(hdr->obs.NChan*sizeof(double));
  RSquared = (double **)malloc(hdr->obs.NChan*sizeof(double));
  ReLconjR = (double **)malloc(hdr->obs.NChan*sizeof(double));
  ImLconjR = (double **)malloc(hdr->obs.NChan*sizeof(double));
  SampleCount = (int **)malloc(hdr->obs.NChan*sizeof(int));
  
  for(i=0;i<hdr->obs.NChan;i++){
    LSquared[i] = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&LSquared[i][0], NPtsProf);
    RSquared[i] = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&RSquared[i][0], NPtsProf);
    ReLconjR[i] = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ReLconjR[i][0], NPtsProf);
    ImLconjR[i] = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ImLconjR[i][0], NPtsProf);
    SampleCount[i] = (int *)malloc(NPtsProf*sizeof(int));
    IZero(&SampleCount[i][0], NPtsProf);
  }


  /*   fits_read_key(Fcal, TDOUBLE, "DUMPMIDSECS",  &DumpMiddleSecs, 
       NULL, &status); status = 0; */
  /*   fits_read_key(Fcal, TDOUBLE, "DUMPREFPER",   &DumpRefPeriod,  
       NULL, &status); status = 0; */
  /*   fits_read_key(Fcal, TDOUBLE, "DUMPREFPHASE", &DumpRefPhase,   
       NULL, &status); status = 0;  */


  colnum = 0;

  if(RunMode->Verbose){
    printf("\n\nCAL FACTORS:\n");
    printf("============\n\n");
    printf(" Freq (MHz)      L^2        R^2         Re(L*R)       Im(L*R) \n");
    printf("-----------     -----      -----       ---------     ---------\n\n");
    fflush(stdout);
  }

  for(i=0;i<hdr->obs.NChan;i++) {

    
    /* Read in data for this channel */
    fits_read_col(Fcal, TDOUBLE, ++colnum, 1, 1, NPtsProf, NULL, 
		  &LSquared[i][0], &anynull, &status); 
    fits_read_col(Fcal, TDOUBLE, ++colnum, 1, 1, NPtsProf, NULL, 
		  &RSquared[i][0], &anynull, &status); 
    fits_read_col(Fcal, TDOUBLE, ++colnum, 1, 1, NPtsProf, NULL, 
		  &ReLconjR[i][0], &anynull, &status); 
    fits_read_col(Fcal, TDOUBLE, ++colnum, 1, 1, NPtsProf, NULL, 
		  &ImLconjR[i][0], &anynull, &status); 
    fits_read_col(Fcal, TLONG,   ++colnum, 1, 1, NPtsProf, NULL, 
		  &SampleCount[i][0], &anynull, &status); 
    
/* Normalize */
    for(j=0;j<NPtsProf;j++){
      LSquared[i][j] /= (double)SampleCount[i][j];
      RSquared[i][j] /= (double)SampleCount[i][j];
      ReLconjR[i][j] /= (double)SampleCount[i][j];
      ImLconjR[i][j] /= (double)SampleCount[i][j];      
    } 
    



    /* Now find cal pulse height for each array, and calculate cal factor in 
       Jy/count for each channel */

    /*  CalFactor[i][0]: corresponds to LSquared[i]
     *	CalFactor[i][1]: corresponds to RSquared[i]
     *	CalFactor[i][2]: corresponds to ReLconjR[i] 
     *	CalFactor[i][3]: corresponds to ImLconjR[i]  */


    printf("GETCAL: Channel %d, %lf MHz:\n",i,hdr->obs.ChanFreq[i]);
    if (RunMode->Swap){

    GetCalFac(hdr, RunMode, LSquared[i], (int)NPtsProf, RunMode->CalValue[1], 
	      &CalFactor[i][0], &Tsys, &Gain);  
    printf("GETCAL: Recvr %s, Gain=%lf K/Jy:\n", hdr->gen.FEName, Gain);
    printf("GETCAL:     LEFT/X: TSYS = %lf Kelvin\n",Tsys);  
    GetCalFac(hdr, RunMode, RSquared[i], (int)NPtsProf, RunMode->CalValue[0], 
	      &CalFactor[i][1], &Tsys, &Gain);
    printf("GETCAL:    RIGHT/Y: TSYS = %lf Kelvin\n",Tsys);  fflush(stdout);
    }
    else {

    GetCalFac(hdr, RunMode, LSquared[i], (int)NPtsProf, RunMode->CalValue[0], 
	      &CalFactor[i][0], &Tsys, &Gain);  
    printf("GETCAL: Recvr %s, Gain=%lf K/Jy:\n", hdr->gen.FEName, Gain);
    printf("GETCAL:     LEFT/X: TSYS = %lf Kelvin\n",Tsys);  fflush(stdout);
    GetCalFac(hdr, RunMode, RSquared[i], (int)NPtsProf, RunMode->CalValue[1], 
	      &CalFactor[i][1], &Tsys, &Gain);
    printf("GETCAL:    RIGHT/Y: TSYS = %lf Kelvin\n",Tsys);  fflush(stdout);

    }

    CalFactor[i][2] = CalFactor[i][3] = sqrt(CalFactor[i][0] * CalFactor[i][1]);

    
    if(RunMode->Verbose)
      printf("    %.2f     %.4f     %.4f     %.4f    %.4f\n",
	     hdr->obs.ChanFreq[i],CalFactor[i][0],CalFactor[i][1],
	     CalFactor[i][2],CalFactor[i][3]);fflush(stdout);



  }

  if(RunMode->Verbose)
    printf("\n");
    
  for(i=0;i<hdr->obs.NChan;i++) {
    free(LSquared[i]); 
    free(RSquared[i]); 
    free(ReLconjR[i]);
    free(ImLconjR[i]);
  }

  free(LSquared);
  free(RSquared); 
  free(ReLconjR);
  free(ImLconjR);
  
    
  


}
