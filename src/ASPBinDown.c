/* A binning down utility for asp-created ascii files */
/* -- RDF 08 Feb 2006 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "BinDownCmdLine.h" 

int main(int argc, char **argv)
{

  int linenum, bin[NBINMAX];
  int NBinsIn;
  float IMJD,MidSecs,Period,Freq,DM,Phase,Wgt;
  float lin[NBINMAX], phi[NBINMAX], phierr[NBINMAX], smptype[NBINMAX];
  char Outfile[256];
  char line[256], Headerline[256];
  char ObsCode[4],junk[32];
  struct RunVars RunMode;
  struct StdProfs InputProfs, OutputProfs;

  FILE *fpin;

  Cmdline *Cmd;

  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /*  if(argc != 2){
    printf("Please only one file as argument and that's it.\n");fflush(stdout);
    exit(1);
    } */

#if 0
  if(Cmd->m4headP && Cmd->NoHeaderP){
    printf("\n Can't choose -m4head AND -nohead!  geez...\n");fflush(stdout);
    exit(1);
  }
  if((Cmd->NoHeaderP && !Cmd->SourceP) || (Cmd->NoHeaderP && !Cmd->FreqP)){
    printf("\n If you choose no header then you must provide source name \n");
    printf(" and observing frequency.  Exiting...\n");fflush(stdout);
    exit(1);
  }
#endif

  strcpy(Outfile,"BinDown.out");

  /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  strcpy(RunMode.Infile,Cmd->Infile); 
  RunMode.Verbose=Cmd->VerboseP;
  RunMode.NBinsOut=Cmd->NBinsOut;

  if((fpin = fopen(Cmd->Infile,"r"))==NULL){
    printf("Could not read file %s. Exiting...\n",Cmd->Infile);fflush(stdout);
    exit(1);
  }

  /* Set rest ofRunMode structure elements to be compatible with ASPFitsReader 
     code */
  RunMode.BinDown=1;
  RunMode.Header=1;
  /* Read in file starting on second line */

  linenum=-1;

#if 0
  if(Cmd->NoHeaderP){ // if there is no header line...
    printf("\nNo header mode chosen.  Will create mark4-style header.\n");
    strcpy(RunMode.Source,Cmd->Source);      
    sprintf(Headerline,"%4.0f. %4d   %7s            %.3f",
	   Cmd->Freq,RunMode.NBinsOut,RunMode.Source,1.0);
    linenum++;
  }
#endif

  while ((fgets(line,128,fpin) != NULL)){
    if(linenum<0){
      if (Cmd->m4headP){ // m4-style header in
	sscanf(line,"%f%d%s%f",&Freq,&NBinsIn,RunMode.Source,&Wgt);
	sprintf(Headerline,"%4.0f. %4d   %7s            %.3f",
		Freq,RunMode.NBinsOut,RunMode.Source,Wgt);
      }
      else{ // asp-pspmtoa style header line
	sscanf(line,"%s%f%f%f%s%f%f%d%s%s%s%f",junk,
	       &IMJD,&MidSecs,&Period,junk,&Freq,&DM,&NBinsIn,
	       ObsCode,junk,RunMode.Source,&Phase);
	sprintf(Headerline,
		"# %.1f %.7f %.10f %ld %.3f %.3f %d %s %d %s %.10f",
		IMJD,MidSecs,Period,(long)1,Freq,DM,RunMode.NBinsOut,
		ObsCode,1,RunMode.Source,Phase);
      }
      linenum++;
    }
    else{
      sscanf(line,"%d%f%f%f%f%f%f%f%f",
	     &bin[linenum],&InputProfs.rstds[linenum],
	     &InputProfs.rstdq[linenum],&InputProfs.rstdu[linenum],
	     &InputProfs.rstdv[linenum],
	     &lin[linenum],&phi[linenum],&phierr[linenum],&smptype[linenum]);
      linenum++;
    }
  }
  fclose(fpin);


  RunMode.NBins=linenum;
  if(RunMode.NBins < RunMode.NBinsOut){
    printf("\nNumber of Bins to output MUST be less than number \n");
    printf("of input bins.  Exiting...\n");fflush(stdout);
    exit(1);
  }
  
  printf("\nInput profile data file: %s\n\n",Cmd->Infile);
  printf("Number of bins in input profile : %d\n\n",RunMode.NBins);
  printf("Number of bins in output profile: %d\n\n",RunMode.NBinsOut);

  



  BinDown(&RunMode,&InputProfs,&OutputProfs);

  //  printf("Header:\n%s\n",Headerline);fflush(stdout);

  MakePol(&RunMode,RunMode.NBins,&OutputProfs);

  WriteStokes(&RunMode,&OutputProfs,Headerline,Outfile);

  exit(0);

}
