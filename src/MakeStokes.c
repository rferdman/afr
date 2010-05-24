/* Calculates Stokes parameters from  L^2, R^2, etc. */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

void MakeStokes(struct ASPHdr *hdr, struct RunVars *RunMode, 
		struct StdProfs *StdProfiles, 
		double *ASquared, double *BSquared,
		double *ReAconjB, double *ImAconjB,
		double *JyPerCount)
{

  int i;
  float TempProf[NBINMAX];
  double Duty,FinalMask[RunMode->NBins];
  double LBase,RBase,RScale,CScale;
  double Lrms,Rrms;

  // float MeanA, MeanB;

  /* Here, scale each array to rms of ASquared */
  if(RunMode->Scale) {
    
    //printf("Scaling each polarization by baseline RMS.\n\n");
    // MeanA=MeanB=0;
    /* flaot-friendly format: */
    for (i=0;i<RunMode->NBins;i++) {
      TempProf[i] = (float)(ASquared[i]);
      // MeanA += (float)ASquared[i];
	
    }    
    Duty = DutyLookup(RunMode->Source);
/*     BMask(&ASquared[0],&RunMode->NBins,&Duty,FinalMask); */
    BMask(&TempProf[0],&RunMode->NBins,&Duty,FinalMask);
/*     Baseline(&ASquared[0],FinalMask,&RunMode->NBins,&LBase,&Lrms); */
    Baseline(&TempProf[0],FinalMask,&RunMode->NBins,&LBase,&Lrms);

    /* float-friendly format: */
    for (i=0;i<RunMode->NBins;i++) {
      TempProf[i] = (float)(BSquared[i]);
      // MeanB += (float)BSquared[i]; 
    }

    Baseline(&TempProf[0],FinalMask,&RunMode->NBins,&RBase,&Rrms);
    RScale = Lrms/Rrms;
    CScale = sqrt(RScale);

    /* Calculate means */
    /*     MeanA /= (float)RunMode->NBins;
	   MeanB /= (float)RunMode->NBins; */

    for(i=0;i<RunMode->NBins;i++) {

      /* Same as other scaling, but each pol separately instead of 
	 relative to L */
      ASquared[i] *= Lrms; /* because it's |L|^2 */
      BSquared[i] *= Rrms; /* because it's |R|^2 */
      ReAconjB[i] *= sqrt(Lrms*Rrms); /* because it's Re(L*R) */
      ImAconjB[i] *= sqrt(Lrms*Rrms); /* because it's Re(L*R) */

#if 0
      BSquared[i] *= RScale;  /* because it's |R|^2 */
      ReAconjB[i] *= CScale;  /* because it's Re(L*R) */
      ImAconjB[i] *= CScale;  /* because it's Im(L*R) */
#endif

    }
    if(RunMode->Verbose)  
      printf("Lrms = %f, Rrms = %f, RScale = %f, Cscale = %f\n",
	     Lrms,Rrms,RScale,CScale);fflush(stdout);
  }

  /* Calculate Stokes Parameters */

  for(i=0;i<RunMode->NBins;i++) {
    /* Calibrate */
    ASquared[i] *= JyPerCount[0];
    BSquared[i] *= JyPerCount[1];
    ReAconjB[i] *= JyPerCount[2];
    ImAconjB[i] *= JyPerCount[3];

    /* if GASP data and lower sideband, do Im(L*R) --> -Im(L*R)*/
    if(!strcmp(hdr->gen.BEName, "xASP"))
       ImAconjB[i] *= RunMode->Sideband;

    /* Make Stokes parameters, depending on whether linear or circular basis */
    StdProfiles->rstds[i] = (float)(ASquared[i]+BSquared[i])/2.;

    if(!strcmp("C",hdr->gen.FEPol)){  /* dual circular */

      if(RunMode->Swap){
/* 	StdProfiles->rstdq[i] = (float)(2.*ReAconjB[i]); */
/* 	StdProfiles->rstdu[i] = (float)(-2.*ImAconjB[i]); */
/* 	StdProfiles->rstdv[i] = (float)(BSquared[i]-ASquared[i]); */
	StdProfiles->rstdq[i] = (float)(ReAconjB[i]);
	StdProfiles->rstdu[i] = (float)(-1.0*ImAconjB[i]);
	StdProfiles->rstdv[i] = (float)(BSquared[i]-ASquared[i])/2.;
      }
      else{
/* 	StdProfiles->rstdq[i] = (float)(2.*ReAconjB[i]); */
/* 	StdProfiles->rstdu[i] = (float)(2.*ImAconjB[i]); */
/* 	StdProfiles->rstdv[i] = (float)(ASquared[i]-BSquared[i]); */
	StdProfiles->rstdq[i] = (float)(ReAconjB[i]);
	StdProfiles->rstdu[i] = (float)(ImAconjB[i]);
	StdProfiles->rstdv[i] = (float)(ASquared[i]-BSquared[i])/2.;
      }
      /*     printf("L,R,L,R = %lf. %lf. %lf. %lf\n",ASquared[i],BSquared[i],
	     ReAconjB[i],ImAconjB[i]);fflush(stdout);  */
      /*     printf("S,Q,U,V = %lf, %lf, %lf, %lf\n",StdProfiles->rstds[i],
	     StdProfiles->rstdq[i],StdProfiles->rstdu[i],
	     StdProfiles->rstdv[i]);fflush(stdout); */
    }
    else if(!strcmp("L",hdr->gen.FEPol)){  /* dual linear */

      if(RunMode->Swap){
/* 	StdProfiles->rstdq[i] = (float)(BSquared[i]-ASquared[i]); */
/* 	StdProfiles->rstdu[i] = (float)(2.*ReAconjB[i]); */
/* 	StdProfiles->rstdv[i] = (float)(2.*ImAconjB[i]); */
	StdProfiles->rstdq[i] = (float)(BSquared[i]-ASquared[i])/2.;
	StdProfiles->rstdu[i] = (float)(ReAconjB[i]);
	StdProfiles->rstdv[i] = (float)(-1.0*ImAconjB[i]);
      }
      else{

/* 	StdProfiles->rstdq[i] = (float)(ASquared[i]-BSquared[i]); */
/* 	StdProfiles->rstdu[i] = (float)(2.*ReAconjB[i]); */
/* 	StdProfiles->rstdv[i] = (float)(-2.*ImAconjB[i]); */
 	StdProfiles->rstdq[i] = (float)(ASquared[i]-BSquared[i])/2.;
	StdProfiles->rstdu[i] = (float)(ReAconjB[i]);
	StdProfiles->rstdv[i] = (float)(ImAconjB[i]);
     }

    }
    else {
      printf("Cannot recognize polarization mode description!  Exiting...\n");
      fflush(stdout);
      exit(1);
      
    }

  }


}
