/* Rotates a profile BY a given phase */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"
#include "RotateCmdLine.h"

int main(int argc, char **argv)
{

  int linenum, bin[NBINMAX];
  float lin[NBINMAX], phi[NBINMAX], phierr[NBINMAX], smptype[NBINMAX];
  struct RunVars RunMode;
  struct StdProfs Profile;
  char line[256], Headerline[256], Outfile[32];
  FILE *fpin;

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
  RunMode.Header=1;

  strcpy(Outfile,"RotateProf.out");
  
  if((fpin = fopen(Cmd->Infile,"r"))==NULL){
    printf("Could not read file %s. Exiting...\n",Cmd->Infile);fflush(stdout);
    exit(1);
  }

  linenum=-1;
  while ((fgets(line,128,fpin) != NULL)){
    if(linenum<0){
      strcpy(Headerline,line);
      linenum++;
   }
    else{
      sscanf(line,"%d%f%f%f%f%f%f%f%f",
	     &bin[linenum],&Profile.rstds[linenum],
	     &Profile.rstdq[linenum],&Profile.rstdu[linenum],
	     &Profile.rstdv[linenum],
	     &lin[linenum],&phi[linenum],&phierr[linenum],&smptype[linenum]);
      linenum++;
    }
  }
  fclose(fpin);

  RunMode.NBins = RunMode.NBinsOut =  linenum;
       
  printf("NBins = %d\n",RunMode.NBins);

  RotateProf(&RunMode, &Profile, Cmd->ByAngle);
    
  MakePol(&RunMode,RunMode.NBins,&Profile);


  WriteStokes(&RunMode,&Profile,Headerline,Outfile);

  exit(0);

}
