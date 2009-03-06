/* Bins down input profile to a profile with less bins*/

#include "ASPCommon.h"

void BinDown(struct RunVars *RunMode, 
	      struct StdProfs *InputProfs,
	      struct StdProfs *OutputProfs)

{

  int    i,j,k;
  int    Faci,Facm;
  double Facf,Point,Fac,Facl,Facr;
  

  FZero(OutputProfs->rstds,NBINMAX);
  FZero(OutputProfs->rstdq,NBINMAX);
  FZero(OutputProfs->rstdu,NBINMAX);
  FZero(OutputProfs->rstdv,NBINMAX);
  /* Here starts the actual binning down part */

  if(RunMode->BinDown) {

    /* Bin-down factor */
    Fac = (double)RunMode->NBins/((double)RunMode->NBinsOut); 
    Faci = (int)Fac;     /* Integer part */
    Facf = Fac - Faci;   /* decimal part */

    if(RunMode->Verbose) printf("Fac, faci, facf: %f %d %f\n",Fac,Faci,Facf);

    Point = 0.;  
    for(i=0;i<RunMode->NBinsOut;i++) {

      j = (int)Point;
      Facm = 0.;
      Facl = 1.0 - (Point - j);
      Facr = Fac - Facl;
      if(Facr>=1.0) Facm = (int)Facr;
      Facr -= Facm;
      /*    printf("i, Facl, Facm, Facr, point: %d %f %d %f %f\n",i, Facl, 
	    Facm, Facr, Point); */
      OutputProfs->rstds[i] += Facl*InputProfs->rstds[j];
      OutputProfs->rstdq[i] += Facl*InputProfs->rstdq[j];
      OutputProfs->rstdu[i] += Facl*InputProfs->rstdu[j];
      OutputProfs->rstdv[i] += Facl*InputProfs->rstdv[j];
      if(i==0){
      }
      for(k=1;k<=Facm;k++) {
	OutputProfs->rstds[i] += InputProfs->rstds[j+k];
	OutputProfs->rstdq[i] += InputProfs->rstdq[j+k];
	OutputProfs->rstdu[i] += InputProfs->rstdu[j+k];
	OutputProfs->rstdv[i] += InputProfs->rstdv[j+k];
      }
      if(i==0){
      }
      Point += Fac;
      j = (int)Point;
      /*    printf("point, j = %f %d\n",point,j); */
      if(j<RunMode->NBins) {
	OutputProfs->rstds[i] += Facr*InputProfs->rstds[j];
	OutputProfs->rstdq[i] += Facr*InputProfs->rstdq[j];
	OutputProfs->rstdu[i] += Facr*InputProfs->rstdu[j];
	OutputProfs->rstdv[i] += Facr*InputProfs->rstdv[j];
      }

      OutputProfs->rstds[i] /= Fac; /* Normalize */
      OutputProfs->rstdq[i] /= Fac;
      OutputProfs->rstdu[i] /= Fac;
      OutputProfs->rstdv[i] /= Fac;
      if(i==0){
      }
    }
 
  } 
  else {  /* No Binning Down */
    for(i=0;i<RunMode->NBinsOut;i++) {
      OutputProfs->rstds[i] = InputProfs->rstds[i]; 
      OutputProfs->rstdq[i] = InputProfs->rstdq[i];
      OutputProfs->rstdu[i] = InputProfs->rstdu[i];
      OutputProfs->rstdv[i] = InputProfs->rstdv[i];
    }
  }

  /* All binned down now */

}
