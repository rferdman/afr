#include "ASPCommon.h"

void WriteStokes(struct RunVars *RunMode,
		 struct StdProfs *Profiles,
		 char *HeadLine,
		 char *Stokesfile)
{

  int i;
  FILE *fpout;
  double x,ptype;

  if ((fpout = fopen(Stokesfile,"w")) == 0)
    { printf("Cannot open %s. Exiting...\n",Stokesfile); exit(1); }


  if(RunMode->Header){ /* Choosing header line */

    fprintf(fpout,"%s\n",HeadLine);
  }
  else{
    fprintf(fpout,"%4.0f. %4d   %7s            1.000\n",
	    RunMode->FSky,RunMode->NBinsOut, RunMode->Source);
  }

  /*  Duty = DutyLookup(RunMode->Source);
  BMask(Profiles->rstds,&RunMode->NBinsOut,&Duty,FinalMask);  
  Baseline(Profiles->rstds,FinalMask,&RunMode->NBinsOut,
  &Savg,&Profiles->Srms); */

  
  //printf("WRITESTOKES:  Srms = %f\n",Profiles->Srms);fflush(stdout);

  for(i=0;i<RunMode->NBinsOut;i++) {
/* see how strong the linear polarization is */
    x = Profiles->stdlin[i]*Profiles->Srms; 
    ptype = 43.1;
    if (x > 1.) ptype=43.2;
    if (x > 2.) ptype=43.3;
    if (x > 3.) ptype=43.4;
    if (x > 4.) ptype=43.5;
    if (x > 5.) ptype=43.6;
    fprintf(fpout,"%5d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%6.1f\n",i,
	    Profiles->rstds[i],Profiles->rstdq[i],
	    Profiles->rstdu[i],Profiles->rstdv[i],
	    Profiles->stdlin[i],Profiles->stdphi[i]*360.0/TWOPI, /* degrees */
	    Profiles->stdphierr[i]*180.0/TWOPI,ptype);
  }
  fclose(fpout);

}

