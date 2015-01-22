/* ========================================================================= */
/* GetPhases                                                                 */
/*     Will take in a par file and some input options, and will output a     */
/*     series of MJDs and specify the spin phase of the pulsar at each of    */
/*     these dates.                                                          */
/*                                                                           */
/* R. Ferdman, 2011 May 23                                                   */
/* ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "polyco.h"
#include "PhaseCmdLine.h"

int FileExists(char *);
int PhaseCalc(struct Polyco *, int, double, double, double *, double *);


int main(int argc, char **argv)
{

  int            n_poly=0, poly_set_used;
  double         mjd, HoursElapsed, RefPhase, RefFreq, RefPeriod;
  char           tempo_cmd[256];
  struct Polyco  Polyco[MAX_PC_SETS];
  Cmdline        *Cmd;

  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);
  /*  showOptionValues();  */

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;
                                                                             
  /* Store program name */
  //strcpy(ProgName,argv[0]);  

  /* Make sure parfile exists */
  if(FileExists(Cmd->ParFile)){
    if (Cmd->VerboseP)  
      printf("Creating polycos for PSR %s from input parameter file %s.\n", 
	     Cmd->PSRName, Cmd->ParFile); 
  }
  else{
    fprintf(stderr,"Could not open file %s.\n", Cmd->ParFile);
    return -1;
  }


  /* If we made it this far, we have a par file we can now use to create 
     a polyco file. */

  /* Open final output file and close right away to create it -- sort of 
     like a "touch" command, in effect */
  /*  if((Fpoly=fopen("poly_final.dat", "w"))==NULL){
    fprintf(stderr,"Could not open file poly_final.dat.  \n");
    return -1;
  }
  fclose(Fpoly); */

  /* Now create polycos for the observing times, timespan of validity, and 
     observing frequency specified by user */

  sprintf(tempo_cmd,
       "tempo -f %s -Zpsr=%s -Zfreq=%.3lf -Zstart=%.3lf -Ztobsh=%.2lf -Zspan=%d -Zsite=%s",
	  Cmd->ParFile, Cmd->PSRName, Cmd->ObsFreq, 
	  /* Start polyco set to be 1 hour early, (and so end 1 hour later) */
	  Cmd->MJDStart-0.04  /* start one hour early */,
	  /* allow an extra 2 hours of polycos */
	  Cmd->TObsHours + 2., /* End two hour later */ 
	  /* default at 30 minutes valid span */
	  Cmd->NSpan, Cmd->Site);
  
  if (Cmd->VerboseP) printf("tempo command:\n%s\n\n\n", tempo_cmd);
  system(tempo_cmd);
  

  /* Polyco file is now made. Now read it in and save into Polyco array*/


  char testchar[16], testfile[32];
  double testfreq, testmjd;
  strcpy(testfile, "polyco.dat");
  strcpy(testchar, Cmd->PSRName);
  testfreq = Cmd->ObsFreq;
  testmjd = Cmd->MJDStart;

  if (Cmd->VerboseP) { 
    printf("Time step = %.2lf minutes = %.4lf hours = %.5lf days\n",
	   Cmd->TimeStep, Cmd->TimeStep/60.0, Cmd->TimeStep/1440.);fflush(stdout);
    printf("MJDStart = %.4lf, TObsHours = %.3lf\n", 
	   Cmd->MJDStart, Cmd->TObsHours);fflush(stdout);
  }

  if((n_poly=GetPoly("polyco.dat", 
		     Cmd->PSRName,
		     &Polyco[0], 
		     Cmd->ObsFreq,
		     Cmd->MJDStart)) < 1) {
    printf("Could not find polycos for all input parameters from \n");
    printf("par file %s. Exiting...\n", Cmd->ParFile);
    exit(1);
  }


  printf("# MJD          Hours    RefPhase    RefFreq     RefPeriod\n");

  /* Run loop over dates and get out PhaseCalc info */
  for (mjd=Cmd->MJDStart; 
       mjd<(Cmd->MJDStart + Cmd->TObsHours/24.0); 
       mjd+=Cmd->TimeStep/1440.) {
    
    HoursElapsed = (mjd - Cmd->MJDStart)*24.0;

  

    /* Run PhaseCalc */  
    if((poly_set_used = PhaseCalc(Polyco, n_poly, floor(mjd), 
				  (mjd-floor(mjd)), &RefPhase, &RefFreq)) < 0){
      printf("Could not find any polyco sets to match MJD. Exiting...\n");
      exit(2);
    }
    
    if(Cmd->MilliSecsP)
      RefPeriod = 1000./RefFreq;
    else
      RefPeriod = 1.0/RefFreq;

    /* Display PhaseCalc output info */    
    printf("%.6lf   %.4lf   %lf    %lf   %.8lf\n", 
	   mjd, HoursElapsed, RefPhase, RefFreq, RefPeriod);

  }




  exit(0);
}


