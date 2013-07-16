#include <stdio.h>
#include <math.h>
// #include "cpgplot.h"
#include "fitsio.h"
#include "ASPCommon.h"

/* Various tools for retrieving and using calibration data */

/* Give phases where cal turns on/off for each polaization, by first scaling
profiles from each polarization so they are the same height */ 

int GetPhases(struct ASPHdr *hdr, struct RunVars *RunMode, 
	      struct CalVars *CalMode, double **ASquared, 
	      double **BSquared, int *OnBin, int *OffBin){


  /*   int            A2Bin[2][2], B2Bin[2][2], ReBin[2][2], ImBin[2][2];   */
  int     i_phase, bin, omit, nomit, chan, nchans;
  double  APlusB[NBINMAX], MedA[NBINMAX], MedB[NBINMAX];
  double  NormA[NBINMAX], NormB[NBINMAX];
  double  CalHeight, LowAvg;//, HighAvg;
  int     NumNearest=100, TempPhaseBin[2];
  int     PhaseBin[2];  
  int     SmallBin, LargeBin;
  int     N1, N2;
  int     junk;
  int     Range1[2], Range2[2];
//  int     OnBin[2], OffBin[2];
  double  Avg1, Avg2;

  DZero(APlusB,NBINMAX);

  nchans=0;
  for (chan=0;chan<hdr->obs.NChan;chan++){

    omit=0;
    for(nomit=0;nomit<CalMode->NOmit;nomit++){
      if(hdr->obs.ChanFreq[chan] == CalMode->FOmit[nomit]){
	omit=1; break;
      }
    }

    /* Median Filter each profile, find average of top and of bottom, then
       normalize each before adding together */

      
    if(!omit){

      MedianFilter(ASquared[chan], MedA, RunMode->NBins, NumNearest);
      MedianFilter(BSquared[chan], MedB, RunMode->NBins, NumNearest);

      
      if(CalMode->ForcePhase < 0.0) { // finding phases from the data
	if(GetCalPhases(MedA, RunMode->NBins, PhaseBin) < 0){
	  printf("Error finding phase bins for individual Prof, chan %d (%lf MHz).\n",
		 chan, hdr->obs.ChanFreq[chan]) ;
	  printf("Exiting...\n");	fflush(stdout);
	  return -3;
	}         
      }
      else{ // forced phase bin - assumes 0.5 duty cycle
	for (i_phase=0; i_phase<2; i_phase++){
	  PhaseBin[i_phase] = (int)(floor((0.5*i_phase + CalMode->ForcePhase)
					      *((double)RunMode->NBins)));
	  //	  printf("PhaseBin = %d\n",
	  //		 PhaseBin[i_phase]);fflush(stdout);
	  /* Correct to be > 0 and <= NBins-1 */
	  if(PhaseBin[i_phase] < 0)
	    PhaseBin[i_phase] += RunMode->NBins;
	  if(PhaseBin[i_phase] >= RunMode->NBins)
	    PhaseBin[i_phase] -= RunMode->NBins;
	}
	

	/* Now just ensure that bin 0 < bin 1*/
	if(PhaseBin[0] > PhaseBin[1]){
	  junk = PhaseBin[0];
	  PhaseBin[0] = PhaseBin[1];
	  PhaseBin[1] = junk;
	}
	
	/* PhaseBin[0] = (int)(floor(CalMode->ForcePhase
				      *((double)RunMode->NBins+0.5)));
           PhaseBin[1] = (int)(floor(0.5+CalMode->ForcePhase
	                              *((double)RunMode->NBins+0.5))); */
      }
      
      /* Cut off the outside 3 bins from each side of cal ranges */
      Range1[0]=PhaseBin[0]+3; Range1[1]=PhaseBin[1]-3;
      Range2[0]=PhaseBin[1]+3; Range2[1]=PhaseBin[0]-3;
      /* Ensure these are > 0 and < NBins */
      for(i_phase=0; i_phase<2; i_phase++){
	if(Range1[i_phase] < 0)
	  Range1[i_phase] += RunMode->NBins;
	if(Range2[i_phase] < 0)
	  Range2[i_phase] += RunMode->NBins;
	if(Range1[i_phase] >= RunMode->NBins)
	  Range1[i_phase] -= RunMode->NBins;
	if(Range2[i_phase] >= RunMode->NBins)
	  Range2[i_phase] -= RunMode->NBins;	
      }

      //      printf("Range1 = [%d, %d]\nRange2 = [%d, %d]\n",
      //	     Range1[0], Range1[1], Range2[0], Range2[1]);fflush(stdout);

/*      CalHeight = GetCalHeight(MedA, RunMode->NBins,PhaseBin,
			       &LowAvg,&HighAvg);  */
      CalHeight = GetCalHeight(MedA, RunMode->NBins, Range1, Range2,
			       &Avg1, &Avg2);

      //      printf("CalHeight = %lf\n", CalHeight);fflush(stdout);exit(0);

      /* Determine which average is low and which is high */
      if(Avg1<Avg2) 
        LowAvg = Avg1;
      else
        LowAvg = Avg2;

      for (bin=0;bin<RunMode->NBins;bin++)
	NormA[bin] = (ASquared[chan][bin] - LowAvg)/fabs(CalHeight);

      /* Now calculate NormB */

      /* PhaseBin already decided in case of forced phase bin, so only get it 
         this time if phase switch bin needs to be calculated from data */
      if(CalMode->ForcePhase < 0.0) { 
	if(GetCalPhases(MedB, RunMode->NBins, PhaseBin) < 0){
	  printf("Error finding phase bins for individual Prof, chan %d (%lf MHz).\n",
		 chan, hdr->obs.ChanFreq[chan]) ;
	  printf("Exiting...\n");
	  fflush(stdout);
	  return -3;
	}    
      }      
     
      Range1[0]=PhaseBin[0]+3; Range1[1]=PhaseBin[1]-3;
      Range2[0]=PhaseBin[1]+3; Range2[1]=PhaseBin[0]-3;
      /* Ensure these are > 0 and < NBins */
      for(i_phase=0; i_phase<2; i_phase++){
	if(Range1[i_phase] < 0)
	  Range1[i_phase] += RunMode->NBins;
	if(Range2[i_phase] < 0)
	  Range2[i_phase] += RunMode->NBins;
	if(Range1[i_phase] >= RunMode->NBins)
	  Range1[i_phase] -= RunMode->NBins;
	if(Range2[i_phase] >= RunMode->NBins)
	  Range2[i_phase] -= RunMode->NBins;	
      }
 
/*      CalHeight = GetCalHeight(MedB, RunMode->NBins,PhaseBin,
			       &LowAvg,&HighAvg); */
      //  printf("FORCEPHASE= %f, TEMPPHASEBIN[0] = %d, TEMPPHASEBIN[1] = %d\n",CalMode->ForcePhase, PhaseBin[0], PhaseBin[1]);fflush(stdout);

      CalHeight = GetCalHeight(MedB, RunMode->NBins,Range1,Range2,
			      &Avg1,&Avg2);

      /* Determine which average is low and which is high */
      if(Avg1<Avg2) 
        LowAvg = Avg1;
      else
        LowAvg = Avg2;

      for (bin=0;bin<RunMode->NBins;bin++)
	NormB[bin] = (BSquared[chan][bin] - LowAvg)/fabs(CalHeight);

      for (bin=0;bin<RunMode->NBins;bin++){
	APlusB[bin] += NormA[bin] + NormB[bin];
      }	
      nchans++;
    }
  }
  for (bin=0;bin<RunMode->NBins;bin++) APlusB[bin] /= ((float)(2*nchans));

  /* Now find phases where cal turns on and off, this time using 
     averaged cal profile */
  /* Again, if -force option was used, this is already done. */
  if(CalMode->ForcePhase < 0.0) { // forced phase bin
    if(GetCalPhases(APlusB, RunMode->NBins, PhaseBin) < 0){
      printf("Error finding phase bins Exiting...\n");
      fflush(stdout);
      return -3;
    }   
  } 
 
 
  if(RunMode->Verbose)
    printf("\nPHASE BINS : %d  and  %d\n",  PhaseBin[0],PhaseBin[1]);

  /* Now figure out which bin ranges are cal on and off -- low = OFF and
     high = ON, for this total power profile */

  SmallBin = IMin(&PhaseBin[0], 2, &junk);
  LargeBin = IMax(&PhaseBin[0], 2, &junk);

  /* Get quick average for each range */

  N1=N2=0;
  Avg1=Avg2=0.;

  /* This should work whether or not one of the phase ranges wraps in phase */
  for(bin=0; bin<RunMode->NBins;bin++) {
    if(bin>SmallBin+3 && bin<LargeBin-3){
      Avg1 += APlusB[bin];
      N1++;
    }
    else if(bin<SmallBin-3 || bin >LargeBin+3){
      Avg2 += APlusB[bin];
      N2++;
    }
  }
  
  Avg1 /= ((double)N1);
  Avg2 /= ((double)N2); 
 
  if(Avg1 < Avg2) {
    OffBin[0]=SmallBin+3; OffBin[1]=LargeBin-3;
    OnBin[0]=LargeBin+3;  OnBin[1]=SmallBin-3;
  }
  else{
    OnBin[0]=SmallBin+3;  OnBin[1]=LargeBin-3;
    OffBin[0]=LargeBin+3; OffBin[1]=SmallBin-3;
  } 
  /* Ensure these are > 0 and < NBins */
  for(i_phase=0; i_phase<2; i_phase++){
    if(OnBin[i_phase] < 0)
      OnBin[i_phase] += RunMode->NBins;
    if(OffBin[i_phase] < 0)
      OffBin[i_phase] += RunMode->NBins;
    if(OnBin[i_phase] >= RunMode->NBins)
      OnBin[i_phase] -= RunMode->NBins;
    if(OffBin[i_phase] >= RunMode->NBins)
      OffBin[i_phase] -= RunMode->NBins;	
  }

  if(RunMode->Verbose){
    printf("ON Bins:   %d --> %d\n",OnBin[0],OnBin[1]);
    printf("OFF Bins:  %d --> %d\n\n",OffBin[0],OffBin[1]);

    //    printf("Avg1 = %lf, N1 = %d\n",Avg1,N1);
    //    printf("Avg2 = %lf, N2 = %d\n\n",Avg2,N2);fflush(stdout);
  }
    
  return 0;

}

/* Find phases at which cal pulse is on/off */
int GetCalPhases(double *Profile, int NPtsProf, int *BinChange)
{

  int    i;
  double Diff;
  int    NumNearest;
  int    BinDiff, BinDiffThresh;
  int    MaxBin1,  MaxBin2;
  double MaxDiff1, MaxDiff2;
  double MedProf[NBINMAX];



  MaxDiff1 = 0.;
  MaxDiff2 = 0.;
  MaxBin1  = 0;
  MaxBin2  = 0;
  BinChange[0] = BinChange[1] = 0;
      
  /* median smooth entire cal profile */

  NumNearest = 100; // MUST BE EVEN!!!

  MedianFilter(Profile, MedProf, NPtsProf, NumNearest);

  /* finished median smoothing */
  
  /* first find largest differences between count values -- 
     these are the bins to note */

  for (i=0;i<NPtsProf;i++) {
    //    if (i==0)
    if (i<2)
      //      Diff = MedProf[0] - MedProf[NPtsProf-1];
      Diff = MedProf[0] - MedProf[NPtsProf-i-1];
    else
      //      Diff = MedProf[i] - MedProf[i-1];
      Diff = MedProf[i] - MedProf[i-2];


    /* Keep a running max diff, one for up and one for down */
    if( (Diff > 0) && (Diff > MaxDiff1) ){
      MaxDiff1 = Diff;
      if(i==0)
	MaxBin1  = NPtsProf-1;
      else
	MaxBin1  = i-1;
    }
    else if ( (Diff < 0) && (Diff < MaxDiff2) ){
      MaxDiff2 = Diff;
      if(i==0)
	MaxBin2  = NPtsProf-1;
      else
	MaxBin2  = i-1;      
    }

#if 0

    if(fabs(Diff) > fabs(MaxDiff1)){ 
      
      MaxDiff2 = MaxDiff1;
      MaxBin2  = MaxBin1;
      MaxDiff1 = Diff;
      //	MaxBin1  = i;
      if(i==0)
	MaxBin1  = NPtsProf-1;
      else
	MaxBin1  = i-1;
      
    }

    else if (fabs(Diff) > fabs(MaxDiff2)){
      MaxDiff2 = Diff;
      MaxBin2  = i;

    }

#endif


  }
  
  if((MaxDiff1 < 0. && MaxDiff2 < 0.) || (MaxDiff1 > 0. && MaxDiff2 > 0.)){
    printf("ERROR: wrong sign on one of MaxDiff1 (Bin %d) or MaxDiff2 (Bin %d).\n",
	   MaxBin1,MaxBin2);
/*    if(QuickPlot(MedProf, NPtsProf) < 0){
      printf("Could not plot cal prof...\n");fflush(stdout);
    } */
    printf("       Exiting...\n");
    fflush(stdout);
    return -1.;
  }


  BinChange[0] = MaxBin1;
  BinChange[1] = MaxBin2;

  BinDiff = BinChange[0] - BinChange[1];
  BinDiffThresh = (int)(0.40 * (float)NPtsProf);

  /* UpBin and DownBin *should* be separated by half the number of bins in 
   * the profile -- so if the difference between UpBin and DownBin is outside 
   * 10% of half the number of bins, then barf */
  if(abs(BinDiff) - BinDiffThresh > BinDiffThresh ) {
    printf("Error: Size of folded pulse is not about\n");
    printf("       half the number of bins!  Exiting...\n");fflush(stdout);
    return -2.;
  }
 
  //free(MedProf);
  //free(Neighbour);

  return 0;
}
 

/* calculate cal height, knowning on/off phases */
/* double GetCalHeight(double *Profile, int NPtsProf, int *BinChange, 
		    double *LowAvg, double *HighAvg) */
double GetCalHeight(double *Profile, int NPtsProf, int *OnBin, int *OffBin, 
		    double *OnAvg, double *OffAvg)
{
  int    i,ion1,ion2,ioff1,ioff2,j;
//  int    LargeBin, SmallBin;
  int    NOn, NOff;
  double CalOn[NBINMAX], CalOff[NBINMAX];
  double AvgOn, AvgOff, RmsOn, RmsOff, SumOn, SumOff, SDevOn, SDevOff;
  double CalHeight;

  NOn = NOff = 0;
  AvgOn = AvgOff = 0.;
  RmsOn = RmsOff = 0.;
  SumOn = SumOff = 0.;
   
#if 0
  if(BinChange[0]>BinChange[1]){
    LargeBin = BinChange[0];
    SmallBin = BinChange[1];
  }
  else{
    LargeBin = BinChange[1];
    SmallBin = BinChange[0];    
  }
#endif    

  //  PulseInside = (DownBin > UpBin);  
  /* this refers to whether cal pulse is entirely inside the profile 
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

  /* CALCULATE RMS AS WELL!!! */  

  if(OnBin[0] < OnBin[1]){   // this means Cal OFF wraps around phase

    for (j=0;j<NPtsProf;j++){

      if(j>OnBin[0] && j<OnBin[1]){  /* Top of pulse */
        CalOn[NOn++]=Profile[j];
      }
      else if(j>OffBin[0] || j<OffBin[1]){
        CalOff[NOff++]=Profile[j];
      }
    
    }

  }
  else{                      // this means Cal ON wraps around phase

    for (j=0;j<NPtsProf;j++){  

      if(j>OffBin[0] && j<OffBin[1]){  /* Top of pulse */
        CalOff[NOff++]=Profile[j];
      }
      else if(j>OnBin[0] || j<OnBin[1]){
        CalOn[NOn++]=Profile[j];
      }
    
    }

  }

  /* Sort each of the CalOn and CalOff arrays */
  //  DSort(NOn,CalOff);
  //  DSort(NOff,CalOff);
  DSort(NOn,CalOn);
  DSort(NOff,CalOff);

  /* Find Avg and RMS of the central 80% of cal profile values for each of 
     CalOn and CalOff */

  
  //  ion=(int)(0.8*(float)NOn);
  //  ioff=(int)(0.8*(float)NOff);
  ion1=(int)(0.1*(float)NOn);
  ioff1=(int)(0.1*(float)NOff);
  ion2=(int)(0.9*(float)NOn);
  ioff2=(int)(0.9*(float)NOff);
  NOn=NOff=0;

  for(i=ion1;i<ion2;i++){
    AvgOn += CalOn[i];
    RmsOn += CalOn[i]*CalOn[i];
    NOn++;
  }
  AvgOn /= ((double)NOn);  
  SDevOn = sqrt( RmsOn/((double)NOn) - (AvgOn*AvgOn) );
  
  //printf("i1 = %d, N1 = %d, Avg1 = %lf, SDev1 = %lf\n",i1,N1,Avg1,SDev1);

  for(i=ioff1;i<ioff2;i++){
    AvgOff += CalOff[i];
    RmsOff += CalOff[i]*CalOff[i];
    NOff++;
  }
  AvgOff /= ((double)NOff);
  SDevOff = sqrt(RmsOff/((double)NOff) - (AvgOff*AvgOff) );

  //printf("i2 = %d, N2 = %d, Avg2 = %lf, SDev2 = %lf\n",i2,N2,Avg2,SDev2);

  /* Now calculate the REAL average, not counting points that are > three sigma
     away from mean for each of CalOn and CalOff arrays */  

  NOn=NOff=0;
  
  if(OnBin[0] < OnBin[1]){   // this means Cal OFF wraps around phase

    for (j=0;j<NPtsProf;j++){

      if(j>OnBin[0] && j<OnBin[1]){  /* Top of pulse */
	//        if (fabs(Profile[j]-AvgOn) < 3.*RmsOn){
        if (fabs(Profile[j]-AvgOn) < 3.*SDevOn){
	  NOn++;
	  SumOn += Profile[j];
        }
      }
      else if(j>OffBin[0] || j<OffBin[1]){
	//	if (fabs(Profile[j]-AvgOff) < 3.*RmsOff){
	if (fabs(Profile[j]-AvgOff) < 3.*SDevOff){
	  NOff++;
	  SumOff += Profile[j];	
        }
      }
    
    }

  }
  else {                     // this means Cal ON wraps around phase

    for (j=0;j<NPtsProf;j++){

      if(j>OffBin[0] && j<OffBin[1]){  /* Top of pulse */
	//        if (fabs(Profile[j]-AvgOff) < 3.*RmsOff){
        if (fabs(Profile[j]-AvgOff) < 3.*SDevOff){
	  NOff++;
	  SumOff += Profile[j];
        }
      }
      else if(j>OnBin[0] || j<OnBin[1]){
	//        if (fabs(Profile[j]-AvgOn) < 3.*RmsOn){
        if (fabs(Profile[j]-AvgOn) < 3.*SDevOn){
	  NOn++;
	  SumOn += Profile[j];	
        }
      }
    
    }

  }

#if 0
  if(Sum1 > Sum2){
    *HighAvg = Sum1/((double)N1);
    *LowAvg  = Sum2/((double)N2);
  }
  else{
    *HighAvg = Sum2/((double)N2);
    *LowAvg  = Sum1/((double)N1);
  }
#endif    
  
  *OnAvg   = SumOn/((double)NOn);
  *OffAvg  = SumOff/((double)NOff);

  CalHeight = *OnAvg - *OffAvg;
  // CalHeight = fabs(*OnAvg - *OffAvg);

  return CalHeight;
}


/* read cal data from fits file with the option of adding all dumps together */
int GetCalData(struct ASPHdr *CalHdr, struct SubHdr *CalSubHdr, 
	       struct RunVars *RunMode, fitsfile *Fcal,
	       double **ASquared, double **BSquared, 
	       double **ReAconjB, double **ImAconjB)
{

  int    i, i_dump;
  int    status=0;
  int    NColumns, NDumps2Use;
  long   NPtsProf=0;
  //char   *HeadLine[NCHMAX];
  
  long    **SampleCount;



  NDumps2Use = 1;

  if(!strcmp(CalHdr->gen.BEName, "xASP")) 
    fits_movnam_hdu(Fcal, BINARY_TBL, "ASPOUT0", 0, &status);   
  
  /* If wish to use all dumps together, set loop to run over all dumps */
  if (RunMode->AddDumps) NDumps2Use = RunMode->NDumps;

  /* Malloc SampleCount here since we don't need it outside of this routine */
  SampleCount = (long **)malloc(CalHdr->obs.NChan*sizeof(long));

  /* Malloc SubHdr to number of dumps to use */
  CalSubHdr = (struct SubHdr *)malloc(NDumps2Use*sizeof(struct SubHdr));

  /* Malloc HeadLine because C is stupid */
  /** for (i=0;i<CalHdr->obs.NChan;i++){
    if( (HeadLine[i] = (char *)malloc(128)) == NULL) {
      printf("HeadLine malloc'ing failed for i = %d\n",i);fflush(stdout);
      return -1;;
    }
  } **/

  for (i_dump=0;i_dump<NDumps2Use;i_dump++){


    if(!strcmp(CalHdr->gen.BEName, "xASP")) {
      fits_get_num_rows(Fcal, &NPtsProf, &status);status=0; /* find NPtsProf */
    }
    else if (!strcmp(CalHdr->gen.FitsType, "PSRFITS")) {
      NPtsProf = CalHdr->redn.RNBinTimeDump;
    }
    else{
      fprintf(stderr, "GetCalData ERROR: Unrecognized file format.\n");
    }




    if(RunMode->NBins != NPtsProf) {
      printf("ERROR: Number of bins in profile (%ld)\n",NPtsProf);
      printf("       does not agree with header value (%d).  Exiting...\n",
	     RunMode->NBins);
      fflush(stdout);
      return -2;
    }

    fits_get_num_cols(Fcal, &NColumns, &status);status=0;

    /* If newer version of FITS file then move back to first HDU */
    if(!strcmp(CalHdr->gen.BEName, "xASP")) 
      if(!strcmp(CalHdr->gen.HdrVer,"Ver1.0.1")) {
	fits_movrel_hdu(Fcal, -1, NULL, &status);
      }

    if(RunMode->Verbose)
      printf("NChan = %d, NPtsProf = %ld\n", CalHdr->obs.NChan,NPtsProf);

    if (ReadData(CalHdr, CalSubHdr, RunMode, Fcal, i_dump,  NPtsProf, 
		 ASquared, BSquared, ReAconjB, ImAconjB, SampleCount) < 0){
      fprintf(stderr, "ASPCal ERROR: Could not read data from ");
      fprintf(stderr, "file %s (Error occured when attempting to read ",
	      RunMode->Infile);
      fprintf(stderr, "dump %d)", i_dump);
      return -1;      
    }
    

  }

  return 0;
}
