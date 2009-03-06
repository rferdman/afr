/* Rotates a profile BY a given phase */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


void RotateProf(struct RunVars *RunMode, struct StdProfs *InputProfs, 
		double ByPhase)
{

  int    i,j,n2;
/*   int    forward=1,real=1; */
  int    back=-1,complex=0;
  double Pha1;


  /* Fourier Transform the input profiles: */


/*   cprofc(InputProfs->rstds,RunMode->NBins,
     InputProfs->stdamp,InputProfs->stdpha);  */
  cprofc(InputProfs->rstds,RunMode->NBins,
	 InputProfs->stdamps,InputProfs->stdphas);
  cprofc(InputProfs->rstdq,RunMode->NBins,
	 InputProfs->stdampq,InputProfs->stdphaq);
  cprofc(InputProfs->rstdu,RunMode->NBins,
	 InputProfs->stdampu,InputProfs->stdphau);
  cprofc(InputProfs->rstdv,RunMode->NBins,
	 InputProfs->stdampv,InputProfs->stdphav);


  /* ROTATE input profile to phase given in argument */

  Pha1 = ByPhase;   

/*    Pha1 = ToPhase;  */
  if(RunMode->Verbose)
    printf("Rotating input profile at Phase %f by %f rad = %f bins = %f phase\n",
	   InputProfs->stdphas[1]/TWOPI, Pha1, Pha1*RunMode->NBins/TWOPI, 
	   Pha1/TWOPI);
    
  for(i=1;i<RunMode->NBins/2;i++) {
/*     InputProfs->stdpha[i] = fmod((InputProfs->stdpha[i]-(i)*Pha1),TWOPI);  */
    InputProfs->stdphas[i] = fmod((InputProfs->stdphas[i]+(i)*Pha1),TWOPI);
    InputProfs->stdphaq[i] = fmod((InputProfs->stdphaq[i]+(i)*Pha1),TWOPI);
    InputProfs->stdphau[i] = fmod((InputProfs->stdphau[i]+(i)*Pha1),TWOPI);
    InputProfs->stdphav[i] = fmod((InputProfs->stdphav[i]+(i)*Pha1),TWOPI);

  }


/*   printf("totwgt: %f\n",*totwgt);   */
  uncprofc(InputProfs->stdamps,InputProfs->stdphas,RunMode->NBins,
	   &InputProfs->stds[0][0]);
  uncprofc(InputProfs->stdampq,InputProfs->stdphaq,RunMode->NBins,
	   &InputProfs->stdq[0][0]);
  uncprofc(InputProfs->stdampu,InputProfs->stdphau,RunMode->NBins,
	   &InputProfs->stdu[0][0]);
  uncprofc(InputProfs->stdampv,InputProfs->stdphav,RunMode->NBins,
	   &InputProfs->stdv[0][0]);

  ffft_(&InputProfs->stds[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdq[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdu[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdv[0][0], &RunMode->NBins, &back, &complex);

  for(i=0;i<RunMode->NBins;i++) 
    for(j=0;j<2;j++) {
      InputProfs->stds[i][j] /= (float)(RunMode->NBins);
      InputProfs->stdq[i][j] /= (float)(RunMode->NBins);
      InputProfs->stdu[i][j] /= (float)(RunMode->NBins);
      InputProfs->stdv[i][j] /= (float)(RunMode->NBins);
    }

  n2 = RunMode->NBins*2;

  /*  ffft_(&InputProfs->stds[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdq[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdu[0][0], &RunMode->NBins, &back, &complex);
  ffft_(&InputProfs->stdv[0][0], &RunMode->NBins, &back, &complex); */

  for(i=0;i<RunMode->NBins;i++) {
    InputProfs->rstds[i]= InputProfs->stds[i][0]/2.;
    InputProfs->rstdq[i]= InputProfs->stdq[i][0]/2.;
    InputProfs->rstdu[i]= InputProfs->stdu[i][0]/2.;
    InputProfs->rstdv[i]= InputProfs->stdv[i][0]/2.;
  }


  /**********************************************************************/


}
