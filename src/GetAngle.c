#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

void GetAngle(struct RunVars *RunMode, struct StdProfs *Profile, 
	      struct StdProfs *StdProfile)
{

  int i;
/*   int forward = 1, real = 1; */
  int back = -1, complex = 0;
  float Shift,eShift,SNR,eSNR,b,errb;
  int ngood;
  double PhaseShift;
  float PhaTmps[NBINMAX/2+1], PhaTmpq[NBINMAX/2+1];
  float PhaTmpu[NBINMAX/2+1], PhaTmpv[NBINMAX/2+1];
  float StdRots[NBINMAX][2], StdRotq[NBINMAX][2];
  float StdRotu[NBINMAX][2], StdRotv[NBINMAX][2];
  float StdLin[NBINMAX], StdPhi[NBINMAX], ProfLin[NBINMAX], ProfPhi[NBINMAX];
  int StdPk,ProfPk;
  float Stds[NBINMAX],StdPeak,ProfPeak;
  double StdLinAvg,StdLinRms,ProfLinAvg,ProfLinRms;
  double Duty, FinalMask[NBINMAX];
  double Wgt,AngTmp, NewAng;
  double CosSum, SinSum;
  struct StdProfs TempProf;

  int HACK = 0;

  Shift = 0.;
  fftfit_(&Profile->rstds[0],&StdProfile->stdamps[1],&StdProfile->stdphas[1],
	  &RunMode->NBins,&Shift,&eShift,&SNR,&eSNR,&b,&errb,&ngood);
  
  if(Shift < 0.0) Shift += RunMode->NBins;
  PhaseShift = Shift/RunMode->NBins*TWOPI;
  if(PhaseShift < 0.0) PhaseShift += TWOPI;


  CosSum = 0.;
  SinSum = 0.;

  FZero(PhaTmps,NBINMAX/2+1);
  FZero(PhaTmpq,NBINMAX/2+1);
  FZero(PhaTmpu,NBINMAX/2+1);
  FZero(PhaTmpv,NBINMAX/2+1);
  FZero(&StdRots[0][0],2*NBINMAX);
  FZero(&StdRotq[0][0],2*NBINMAX);
  FZero(&StdRotu[0][0],2*NBINMAX);
  FZero(&StdRotv[0][0],2*NBINMAX);
  FZero(StdLin,NBINMAX);
  FZero(StdPhi,NBINMAX);
  FZero(ProfLin,NBINMAX);
  FZero(ProfPhi,NBINMAX);
  
  for(i=1;i<RunMode->NBins/2+1;i++) {
    PhaTmps[i] = fmod((StdProfile->stdphas[i]+(i)*PhaseShift),TWOPI);
    PhaTmpq[i] = fmod((StdProfile->stdphaq[i]+(i)*PhaseShift),TWOPI);
    PhaTmpu[i] = fmod((StdProfile->stdphau[i]+(i)*PhaseShift),TWOPI);
    PhaTmpv[i] = fmod((StdProfile->stdphav[i]+(i)*PhaseShift),TWOPI);
  }

    
  uncprofc(StdProfile->stdamps,PhaTmps,RunMode->NBins,&StdRots[0][0]);
  uncprofc(StdProfile->stdampq,PhaTmpq,RunMode->NBins,&StdRotq[0][0]);
  uncprofc(StdProfile->stdampu,PhaTmpu,RunMode->NBins,&StdRotu[0][0]);
  uncprofc(StdProfile->stdampv,PhaTmpv,RunMode->NBins,&StdRotv[0][0]);
  
  ffft_(&StdRots[0][0],&RunMode->NBins,&back,&complex);
  ffft_(&StdRotq[0][0],&RunMode->NBins,&back,&complex);
  ffft_(&StdRotu[0][0],&RunMode->NBins,&back,&complex);
  ffft_(&StdRotv[0][0],&RunMode->NBins,&back,&complex);
    
  for(i=0;i<RunMode->NBins;i++) {
    Stds[i] = StdRots[i][0];
  }
  FindPeak(Stds, &RunMode->NBins,&StdPk);
  FindPeak(&Profile->rstds[0],&RunMode->NBins,&ProfPk);
  if(RunMode->Verbose)
    printf("Standard peak is at: %d, prof %d\n",StdPk,ProfPk);



  /* Calculate cossum and sinsum */

  /*   for(m=0;m<hdr->obs.NChan;m++){ */
  for(i=0;i<RunMode->NBins;i++){
    StdLin[i] = sqrt(StdRotq[i][0]*StdRotq[i][0] + StdRotu[i][0]*StdRotu[i][0]);
    StdPhi[i] = atan2(StdRotu[i][0],StdRotq[i][0]);
    ProfLin[i] = sqrt(Profile->rstdq[i]*Profile->rstdq[i] +
		      Profile->rstdu[i]*Profile->rstdu[i]);
    ProfPhi[i] = atan2(Profile->rstdu[i],Profile->rstdq[i]);
  }

  Duty = DutyLookup(RunMode->Source);

  BMask(StdProfile->rstds,&RunMode->NBins,&Duty,FinalMask);
  Baseline(StdLin,FinalMask,&RunMode->NBins,&StdLinAvg,&StdLinRms);
  StdPeak = FindPeak(StdLin, &RunMode->NBins, &StdPk);

  BMask(Profile->rstds,&RunMode->NBins,&Duty,FinalMask);
  Baseline(ProfLin,FinalMask,&RunMode->NBins,&ProfLinAvg,&ProfLinRms);
  ProfPeak = FindPeak(ProfLin, &RunMode->NBins,&ProfPk);

  for(i=0;i<RunMode->NBins;i++) {
    StdLin[i] = (StdLin[i]-StdLinAvg)/(StdPeak-StdLinAvg);
    ProfLin[i] = (ProfLin[i]-ProfLinAvg)/(ProfPeak-ProfLinAvg);
  }


  /* Now calculate sines and cosines */
  if(HACK){

    printf("It's A HACK!!!\n");fflush(stdout);

    for(i=30;i<55;i++) {
      Wgt = StdLin[i]*StdLin[i]*ProfLin[i]*ProfLin[i];
      AngTmp = ProfPhi[i] - StdPhi[i];   
      CosSum += Wgt*cos(AngTmp);
      SinSum += Wgt*sin(AngTmp);
    }


  }
  else {

    for(i=0;i<RunMode->NBins;i++) {
      Wgt = StdLin[i]*StdLin[i]*ProfLin[i]*ProfLin[i];
      AngTmp = ProfPhi[i] - StdPhi[i];  
      CosSum += Wgt*cos(AngTmp);
      SinSum += Wgt*sin(AngTmp);
    }

  }


  /* Re-rotate profs here */

  NewAng = -atan2(SinSum,CosSum);  

  if(RunMode->Verbose)
     printf("NewAng = %lf\n",NewAng);fflush(stdout);

  for(i=0;i<RunMode->NBins;i++){
    TempProf.rstds[i] = Profile->rstds[i];
    TempProf.rstdu[i] = (cos(NewAng)*Profile->rstdu[i] + 
 			 sin(NewAng)*Profile->rstdq[i]); 
    TempProf.rstdq[i] = (cos(NewAng)*Profile->rstdq[i] - 
 			 sin(NewAng)*Profile->rstdu[i]); 
      
    TempProf.rstdv[i] = Profile->rstdv[i];

  }
    
  memcpy(Profile,&TempProf,sizeof(struct StdProfs));
    
    



}
