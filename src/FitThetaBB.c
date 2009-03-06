#include <math.h>
#include <stdio.h>
#include <string.h>
// #include "cpgplot.h"
#include "ASPCommon.h"

int GetThetaBB(struct RunVars *RunMode, struct ASPHdr *hdr)
{
  int       i,j;
  int       Diagnose=1;
  int       ThetaFreqs, ChanMatchCount;
  double    Freq[NCHMAX], ThetaBB[NCHMAX];
  char      line[128];


  FILE      *Fin;


  /* Read in ThetaBB file */
  if((Fin = fopen(RunMode->ThetaBBfile,"r")) == NULL){
    printf("Could not open file %s.  Exiting...\n",RunMode->ThetaBBfile);
    fflush(stdout);
    return -1;
  }
  
  i = 0;

  while ((fgets(line,128,Fin) != NULL)){
    sscanf(line,"%lf%lf", &Freq[i], &ThetaBB[i]);
    i++;
  }
  fclose(Fin);
  ThetaFreqs = i;


  /* Now match ThetaBB freq's to pulsar data file freq's... */
  ChanMatchCount = 0;
  for (i=0;i<hdr->obs.NChan;i++){
    
    for (j=0;j<ThetaFreqs;j++){
      if (Freq[j] == hdr->obs.ChanFreq[i]){
	RunMode->ThetaBB[i] = ThetaBB[j];
	ChanMatchCount++;
	break;
      }      
    }
  }
  if (ChanMatchCount != hdr->obs.NChan){
    printf("ThetaBB file %s does not cover all frequencies in pulsar data.\n",
	   RunMode->ThetaBBfile);
    return -2;
  }

  if(Diagnose){
    printf("\nTHETABB input:\n\n");
    for (i=0;i<hdr->obs.NChan;i++)
      printf("%6.1lf --> %lf\n",hdr->obs.ChanFreq[i],RunMode->ThetaBB[i]);
  }

  return 0;

}




int FitThetaBB(struct RunVars *RunMode, struct ASPHdr *hdr, 
	       struct StdProfs *Profile, int chan, int nscan)
{

  int             i;
  double          Wgt;
  struct StdProfs TempProf;

  int Diagnose=1;
  int PeakBin;
  double IPeak;
  FILE  *fptest;
  

  /* Remove Baseline first */

  RemoveBase(RunMode, RunMode->NBins, Profile);
 

  /* Assume equal weight for now */

  Wgt = 1.0;


  /* Apply ThetaBB */

  if(!strcmp("C",hdr->gen.FEPol)){

    for(i=0;i<RunMode->NBins;i++){
      TempProf.rstds[i] = Wgt*Profile->rstds[i];
      TempProf.rstdq[i] = Wgt*(cos(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdq[i]
			       - sin(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdu[i]);
      TempProf.rstdu[i] = Wgt*(cos(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdu[i] 
			       + sin(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdq[i]);
      TempProf.rstdv[i] = Wgt*Profile->rstdv[i];
    }

  }
  else if(!strcmp("L",hdr->gen.FEPol)){

   for(i=0;i<RunMode->NBins;i++){
      TempProf.rstds[i] = Wgt*Profile->rstds[i];
      TempProf.rstdu[i] = Wgt*(cos(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdu[i] 
			       - sin(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdv[i]);
      TempProf.rstdv[i] = Wgt*(cos(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdv[i]
			       + sin(-(RunMode->ThetaBB[chan]))*
			       Profile->rstdu[i]);
      TempProf.rstdq[i] = Wgt*Profile->rstdq[i];
    }

  }
  
  /* Output ThetaBB-corrected file */

  /* If desired, check that ThetaBB is being applied properly */
  if(Diagnose && (nscan == 0)){
    if((fptest = fopen("testUV_thetabb.dat","a")) == NULL){
      printf("Cannot open file testUV_thetabb.dat.  Sorry dude.\n");
      return -1;
    }
    IPeak = FindPeak(Profile->rstds,&RunMode->NBins,&PeakBin);

    fprintf(fptest,"%6.1lf MHz:\n",hdr->obs.ChanFreq[chan]);fflush(fptest);
    fprintf(fptest,"   U:  %6.3f  --- %5.2lf deg --->  %6.3f\n",
	   Profile->rstdu[PeakBin],-(360./TWOPI)*RunMode->ThetaBB[chan],
	    TempProf.rstdu[PeakBin]);fflush(fptest);
    fprintf(fptest,"   V:  %6.3f  --- %5.2lf deg --->  %6.3f\n\n",
	   Profile->rstdv[PeakBin],-(360./TWOPI)*RunMode->ThetaBB[chan],
	    TempProf.rstdv[PeakBin]);fflush(fptest);
    fclose(fptest);
    
  }

  memcpy(Profile,&TempProf,sizeof(struct StdProfs));


  return 0;
}




