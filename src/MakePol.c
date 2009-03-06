/* Bins down input profile to a profile with less bins*/

#include "ASPCommon.h"

void RemoveBase(struct RunVars *RunMode, int NBins, 
		struct StdProfs *InputProfs)

{

  int    i;
  double SBase,VBase,UBase,QBase;
  //  double Srms,Vrms,Linrms,Urms,Qrms;
  double Duty,SPeak,p0218;
  float  s0218[NBINMAX];
  double FinalMask[NBINMAX];
  int    spk;
  static int first=1;
  

  FZero(InputProfs->stdlin,NBINMAX);
  FZero(InputProfs->stdphi,NBINMAX);
  FZero(InputProfs->stdphierr,NBINMAX);

  /* Now do baseline subtraction */

  //  printf("SOURCE = %s\n",RunMode->Source);fflush(stdout);
  Duty = DutyLookup(RunMode->Source);
  BMask(InputProfs->rstds,&NBins,&Duty,FinalMask);
  Baseline(InputProfs->rstds,FinalMask,&NBins,&SBase,&InputProfs->Srms);

  if(RunMode->Verbose) printf("Srms = %f\n",InputProfs->Srms);fflush(stdout);


  Baseline(InputProfs->rstdv,FinalMask,&NBins,&VBase,&InputProfs->Vrms);
  Baseline(InputProfs->rstdq,FinalMask,&NBins,&QBase,&InputProfs->Qrms);
  Baseline(InputProfs->rstdu,FinalMask,&NBins,&UBase,&InputProfs->Urms);

  /* special provision */
  if(strncmp(RunMode->Source,"0218",4) == 0) {  
    printf("Changing 0218 baseline.\n");
    for(i=0;i<NBins;i++) 
      s0218[i] = -InputProfs->rstds[i];
    p0218 = FindPeak(s0218,&NBins,&spk);
    SBase = -p0218;
  }

  /*  printf("U = %lf, UBase = %lf, UBase/U = %lf\n",InputProfs->rstdu[1910],UBase,UBase/InputProfs->rstdu[1910]);fflush(stdout);  
  printf("V = %lf, VBase = %lf, VBase/V = %lf\n",InputProfs->rstdv[1910],VBase,VBase/InputProfs->rstdv[1910]);fflush(stdout);
  printf("V/U = %lf\n\n",InputProfs->rstdv[1910]/InputProfs->rstdu[1910]);fflush(stdout);*/


  if(!RunMode->NoBase){
    for(i=0;i<NBins;i++) {
      InputProfs->rstds[i] -= SBase;
      InputProfs->rstdq[i] -= QBase;
      InputProfs->rstdu[i] -= UBase;
      InputProfs->rstdv[i] -= VBase;  
    }
  }
  else{
    if (first){
      printf("Baseline subtraction turned OFF.\n\n");fflush(stdout);
      first=0;
    }
  }


  /*  printf("U = %lf, UBase = %lf, UBase/U = %lf\n",InputProfs->rstdu[1910],UBase,UBase/InputProfs->rstdu[1910]);fflush(stdout);  
  printf("V = %lf, VBase = %lf, VBase/V = %lf\n",InputProfs->rstdv[1910],VBase,VBase/InputProfs->rstdv[1910]);fflush(stdout);  
  printf("V/U = %lf\n\n",InputProfs->rstdv[1910]/InputProfs->rstdu[1910]);fflush(stdout); */

  SPeak =  FindPeak(InputProfs->rstds,&NBins, &spk);
  InputProfs->SNR = SPeak*InputProfs->Srms;

  


}


void MakePol(struct RunVars *RunMode, int NBins, 
	     struct StdProfs *InputProfs)

{

  int    i;
  double LinAvg;
  //  double SBase,VBase,LinAvg,UBase,QBase;
  //  double Srms,Vrms,Linrms,Urms,Qrms;
  //  double Duty,SPeak,p0218;
  //  float  s0218[NBINMAX];
  double Duty;
  double FinalMask[NBINMAX];
  //  int    spk;
 

  /* Remove baseline */

  RemoveBase(RunMode, NBins, InputProfs);


  /* Calculate linear polarization, position angle, and PA error */
 
  for(i=0;i<NBins;i++) {
    /* Q^2 + U^2 */
    InputProfs->stdlin[i] = sqrt(InputProfs->rstdq[i]*InputProfs->rstdq[i] + 
				  InputProfs->rstdu[i]*InputProfs->rstdu[i]);  
    /* arctan(U/Q) */
    InputProfs->stdphi[i] = -atan2(InputProfs->rstdu[i],InputProfs->rstdq[i])/2.;  
    if(RunMode->FlipPA) InputProfs->stdphi[i] = -InputProfs->stdphi[i];
    InputProfs->stdphierr[i] = 
      sqrt((InputProfs->rstdq[i]*InputProfs->rstdq[i]/
	    (InputProfs->Urms*InputProfs->Urms) + 
	    InputProfs->rstdu[i]*InputProfs->rstdu[i]/
	    (InputProfs->Qrms*InputProfs->Qrms))/
	   (InputProfs->stdlin[i]*InputProfs->stdlin[i]*InputProfs->stdlin[i]*
	    InputProfs->stdlin[i]) );  
  }

  Duty = DutyLookup(RunMode->Source);
  BMask(InputProfs->rstds,&NBins,&Duty,FinalMask);
  Baseline(InputProfs->stdlin,FinalMask,&NBins,&LinAvg,&InputProfs->Linrms);
  for(i=0;i<NBins;i++) {
    InputProfs->stdlin[i] -= LinAvg;
  }

  if(RunMode->Verbose){
    printf("Q, U, lin, total rms: %f %f %f %f\n",
	   1./InputProfs->Qrms, 1./InputProfs->Urms, 
	   1./InputProfs->Linrms,1./InputProfs->Srms);
    printf("SNR: %f\n",InputProfs->SNR);
  }

}           
