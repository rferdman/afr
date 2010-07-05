/* Takes in two profiles, and rotates the second to match the first, based */
/* on shift found by fftfit  */
/*   -- RDF, 08 Feb 2006 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"
#include "MatchCmdLine.h"

int main(int argc, char **argv)
{
  int i, nprof, bin[NBINMAX], NBins;
  int ngood;
  float  profs[NBINMAX],amps[NBINMAX], phas[NBINMAX];
  float Shift,EShift,SNR,ESNR,b,errb;
  double ByAngle;
  char Headerline[256], Outfile[32], RefOutfile[32], Difffile[32];
  struct RunVars RunMode;
  struct StdProfs Profile[2], DiffProf;

  Cmdline *Cmd;

  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  RunMode.Verbose=Cmd->VerboseP;
  RunMode.NoBase=Cmd->NoBaseP;
  RunMode.Header=1;

  strcpy(Outfile,"MatchedProf.out");
  strcpy(Difffile,"DiffProf.out");
  strcpy(RefOutfile,"RefMatchedProf.out");

  /* Profile 0 = standard prof and profile 1 = "to be fitted" prof */

  /* read in each profile in turn */
  for(nprof=0;nprof<2;nprof++){

    if ( ReadASPAsc(Cmd->Infile[nprof], &Headerline[0], bin,  
		    &Profile[nprof], &NBins) < 0) {
      printf("Error in reading file %s.\n",Cmd->Infile[nprof]);
      fflush(stdout);
      exit(1);
    }
       
    if(nprof==0){
      RunMode.NBins = RunMode.NBinsOut = NBins;
      printf("NBins = %d\n",RunMode.NBins);
    }
    else{
      if (NBins != RunMode.NBins){
	printf("\nBoth profiles MUST have same number of bins!  (%d != %d) Exiting...\n",RunMode.NBins,NBins);
	fflush(stdout);
	exit(1);
      }
    }
  }
    
  /* Now have both profiles read in.  cprofc the first of the two profiles */
  cprofc(Profile[0].rstds,RunMode.NBins,
	 Profile[0].stdamps,Profile[0].stdphas);


  memcpy(amps,Profile[0].stdamps,sizeof(float)*NBINMAX);
  memcpy(phas,Profile[0].stdphas,sizeof(float)*NBINMAX);
  memcpy(profs,Profile[1].rstds,sizeof(float)*NBINMAX);

  /* Now fftfit to find shift required in second profile */
  fftfit_(profs,&amps[1],&phas[1],
	  &RunMode.NBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);  

  /* Now shift second profile -- convert Shift from bins to radians */
  ByAngle = (double)(-Shift/RunMode.NBins*TWOPI);

  printf("Rotating %s by %lf radians\n",Cmd->Infile[1],ByAngle);
  printf("Shift: %lf radians, eShift: %lf radians\n",
	 ByAngle,EShift/RunMode.NBins*TWOPI);
  printf("Scale factor: %f +- %f \n",b,errb);
  //  printf("Headerline = %s\n", Headerline);
  
  RotateProf(&RunMode, &Profile[1], ByAngle);

  if (Cmd->ScaleP) {
    //RotateProf(&RunMode, &Profile[0], 0.);

    for(i=0;i<RunMode.NBins;i++){
      Profile[1].rstds[i] /= b;
      Profile[1].rstdq[i] /= b;
      Profile[1].rstdu[i] /= b;
      Profile[1].rstdv[i] /= b;
    } 

    //    MakePol(&RunMode,RunMode.NBins,&Profile[0]);
    //WriteStokes(&RunMode,&Profile[0],Headerline,RefOutfile);
  }

  MakePol(&RunMode,RunMode.NBins,&Profile[1]);

  WriteStokes(&RunMode,&Profile[1],Headerline,Outfile);

  /* Calculate and output difference profile is asked for */
  if(Cmd->DiffP) {
    MakePol(&RunMode,RunMode.NBins,&Profile[0]);
    for(i=0;i<RunMode.NBins;i++){
      DiffProf.rstds[i] = Profile[1].rstds[i] - Profile[0].rstds[i];
      DiffProf.rstdq[i] = Profile[1].rstdq[i] - Profile[0].rstdq[i];
      DiffProf.rstdu[i] = Profile[1].rstdu[i] - Profile[0].rstdu[i];
      DiffProf.rstdv[i] = Profile[1].rstdv[i] - Profile[0].rstdv[i];
    } 
    WriteStokes(&RunMode,&DiffProf,Headerline,Difffile);

  }

  exit(0);

}






