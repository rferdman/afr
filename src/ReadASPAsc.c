/**  function to read in ascii profile file  **/
/**                    -- RDF, 21 July 2008  **/

#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

int ReadASPAsc(char *Infile, char *Headerline, int *bin, 
	       struct StdProfs *Profile, int *NBins) 
{
  
  int linenum, retval=0; 
  float  lin[NBINMAX], phi[NBINMAX], phierr[NBINMAX], smptype[NBINMAX];  
  char line[512], *ptr;
  FILE *fpin;
  
  // printf("Opening file %s\n",Infile);
  
  if((fpin = fopen(Infile,"r"))==NULL){
    printf("Could not read file %s. Exiting...\n",Infile);
    fflush(stdout);
    return -1;;
  }
  linenum=-1;
  while ((fgets(line,256,fpin) != NULL)){
    if(linenum<0){
      //      strcpy(Healine,Headerline);
      strcpy(Headerline,line);
      /* Strip off the newline character */
      if( (ptr = strchr(Headerline, '\n')) != NULL)
	*ptr = '\0';
      linenum++;
    }
    else{
      sscanf(line,"%d%f%f%f%f%f%f%f%f",
	     &bin[linenum],
	     &Profile->rstds[linenum],&Profile->rstdq[linenum],
	     &Profile->rstdu[linenum],&Profile->rstdv[linenum],
	     &lin[linenum],&phi[linenum],&phierr[linenum],&smptype[linenum]);
      linenum++;
    }
  }
  fclose(fpin);
  
  *NBins=linenum;

  // printf("Headerline = %s\n", Headerline);
  
  retval = 0;
  
  return retval;
  
}
