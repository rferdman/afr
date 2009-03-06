#include <math.h>
#include <stdio.h>
#include <string.h>
// #include "cpgplot.h"
#include "ASPCommon.h"

int GetMueller(char *Muellerfile, float *MM, struct ASPHdr *hdr)
{
  int     i,k,chan;
  int     ChanMatch, MatchedChan[NCHMAX];
  int     CheckMueller=1;
  double  Freq;
  char    line[128];
  FILE    *Fin, *fpcheck;


  /* Read in Mueller Matrix file -- NB: Matrix is iverted already! */
  if((Fin = fopen(Muellerfile,"r")) == NULL){
    printf("Could not open file %s.  Exiting...\n",Muellerfile);
    fflush(stdout);
    return -1;
  }
  /* 

  Input Mueller matrix input file format:
     
  FREQ[0]
  MM[0][0]....MM[0][3]
  MM[1][0]....MM[1][3]
  MM[2][0]....MM[2][3]
  MM[3][0]....MM[3][3]
  FREQ[0]
  MM[0][0]....MM[0][3]
  MM[1][0]....MM[1][3]
  MM[2][0]....MM[2][3]
  MM[3][0]....MM[3][3]  
  ...etc...

  */


  //  memset(&MM[0],0,sizeof(float)*NCHMAX*16);

  /* set diagonal elements to 1 */
  for (chan=0;chan<hdr->obs.NChan;chan++){
    MM[16*chan+0] = MM[16*chan+5] = MM[16*chan+10] = MM[16*chan+15] = 1.0;
  }

  ChanMatch=0;

  while ((fgets(line,128,Fin) != NULL)){

    sscanf(line,"%lf", &Freq); 
    for (chan=0;chan<hdr->obs.NChan;chan++){
      if(Freq==hdr->obs.ChanFreq[chan]){

	for(i=0;i<4;i++){
	  fgets(line,128,Fin);
	  sscanf(line,"%f%f%f%f", &MM[16*chan+4*i], &MM[16*chan+4*i+1], 
		 &MM[16*chan+4*i+2], &MM[16*chan+4*i+3] );
	}
	MatchedChan[ChanMatch++] = chan;
	printf("Found a Mueller Matrix for %6.1lf MHz...\n",
	       hdr->obs.ChanFreq[chan]);fflush(stdout);
	break;
      }
     
    }
    /*    sscanf(line,"%f%f%f%f", MM[4*i], MM[4*i+1], 
	  MM[4*i+2], MM[4*i+3] );
	  i++; */


  }
  fclose(Fin);

  if(ChanMatch == 0){
    printf("No frequency channels match those in Mueller matrix file!\n");
    return -1;
  }
  
  if(ChanMatch < hdr->obs.NChan){
    printf("\nWARNING: Mueller matrix NOT found for every channel.\n");
    /* check to see which frequencies were missed in Mueller Matrix file */
    for(k=0;k<ChanMatch;k++){
      for(chan=0;chan<hdr->obs.NChan;chan++){
	if(MatchedChan[k]==chan) break;
	if(chan==hdr->obs.NChan-1)
	  printf("   %lf MHz does NOT have a Mueller matrix\n",
		 hdr->obs.ChanFreq[chan]);	  
      }      
    }
  }
  
  
  /*  if (i!=4) {  
      printf("Your Mueller Matrix file does not have a 4x4 matrix!\n");
      fflush(stdout);
      return -1;
      } */


  /* check that MM read in is correct */
  if (CheckMueller){
    if((fpcheck = fopen("checkmueller.dat","w")) == NULL){
      printf("Error in opening checkmueller.dat...\n");
      return -2;
    }
    for(chan=0;chan<hdr->obs.NChan;chan++){
      fprintf(fpcheck,"%lf\n",hdr->obs.ChanFreq[chan]);
      for(i=0;i<4;i++)
	fprintf(fpcheck,"%9.6f  %9.6f  %9.6f  %9.6f\n",
		MM[16*chan+4*i],   MM[16*chan+4*i+1], 
		MM[16*chan+4*i+2], MM[16*chan+4*i+3]);
    }
    fclose(fpcheck);
  }


  return 0;

}






int FitMueller(struct RunVars *RunMode, struct ASPHdr *hdr, 
	       struct StdProfs *Profile, int chan)
{

  int             i;  
  //  double          MM[16];
  struct StdProfs TempProf;

  int Diagnose=0;
  static int dump=0;
  double SPeak;
  int PeakBin;
  FILE *fptest;

#if 0  
  if (Diagnose){
    /* open file for diagnostic tests */
    if ((fptest = fopen("check_mueller.dat","w")) == NULL) {
      printf("Cannot open check_mueller.dat. Exiting...\n"); 
      exit(1);
    }

    /* print mueller matrix */
    if(dump==1){
      printf("%7.4f %7.4f %7.4f %7.4f\n",
	     RunMode->MM[0], RunMode->MM[1], RunMode->MM[2], RunMode->MM[3]);
      printf("%7.4f %7.4f %7.4f %7.4f\n",
	     RunMode->MM[4], RunMode->MM[5], RunMode->MM[6], RunMode->MM[7]);
      printf("%7.4f %7.4f %7.4f %7.4f\n",
	     RunMode->MM[8], RunMode->MM[9], RunMode->MM[10],RunMode->MM[11]);
      printf("%7.4f %7.4f %7.4f %7.4f\n",
	     RunMode->MM[12],RunMode->MM[13],RunMode->MM[14],RunMode->MM[15]);
      fflush(stdout);
    }
    fclose(fptest);

    if ((fptest = fopen("check_muel_stokes.dat","w")) == 0) {
      printf("Cannot open check_muel_stokes.dat. Exiting...\n"); 
      exit(1);
    }
  }
#endif


  //  char           ProgName[256];
  


  /* Remove Baseline first */

  //  RemoveBase(RunMode, RunMode->NBins, Profile);
 
  /* Zero out TempProf array */
  memset(&TempProf,0,sizeof(struct StdProfs));

  /* Diagnostic stuff -- find peak, print out MM, and find 
   * effect of each element on matrix multiplication */    

  if(Diagnose)
    SPeak =  FindPeak(Profile->rstds, &RunMode->NBins, &PeakBin);

  /* Apply Mueller Matrix */

  if(!strcmp("C",hdr->gen.FEPol)){
 
    for(i=0;i<RunMode->NBins;i++){
      TempProf.rstds[i] = (RunMode->MM[16*chan + 0]*Profile->rstds[i] 
			   + RunMode->MM[16*chan + 1]*Profile->rstdv[i]
			   + RunMode->MM[16*chan + 2]*Profile->rstdu[i]
			   - RunMode->MM[16*chan + 3]*Profile->rstdq[i]);
      TempProf.rstdv[i] =  (RunMode->MM[16*chan + 4]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 5]*Profile->rstdv[i]
			    + RunMode->MM[16*chan + 6]*Profile->rstdu[i]
			    - RunMode->MM[16*chan + 7]*Profile->rstdq[i]);
      TempProf.rstdu[i] =  (RunMode->MM[16*chan + 8]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 9]*Profile->rstdv[i]
			    + RunMode->MM[16*chan + 10]*Profile->rstdu[i]
			    - RunMode->MM[16*chan + 11]*Profile->rstdq[i]);
      TempProf.rstdq[i] =  (RunMode->MM[16*chan + 12]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 13]*Profile->rstdv[i]
			    + RunMode->MM[16*chan + 14]*Profile->rstdu[i]
			    - RunMode->MM[16*chan + 15]*Profile->rstdq[i]);
    }
    
  }
  else if(!strcmp("L",hdr->gen.FEPol)){
    for(i=0;i<RunMode->NBins;i++){
      TempProf.rstds[i] = (RunMode->MM[16*chan + 0]*Profile->rstds[i] 
			   + RunMode->MM[16*chan + 1]*Profile->rstdq[i]
			   + RunMode->MM[16*chan + 2]*Profile->rstdu[i]
			   + RunMode->MM[16*chan + 3]*Profile->rstdv[i]);
      TempProf.rstdq[i] =  (RunMode->MM[16*chan + 4]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 5]*Profile->rstdq[i]
			    + RunMode->MM[16*chan + 6]*Profile->rstdu[i]
			    + RunMode->MM[16*chan + 7]*Profile->rstdv[i]);
      TempProf.rstdu[i] =  (RunMode->MM[16*chan + 8]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 9]*Profile->rstdq[i]
			    + RunMode->MM[16*chan + 10]*Profile->rstdu[i]
			    + RunMode->MM[16*chan + 11]*Profile->rstdv[i]);
      TempProf.rstdv[i] =  (RunMode->MM[16*chan + 12]*Profile->rstds[i] 
			    + RunMode->MM[16*chan + 13]*Profile->rstdq[i]
			    + RunMode->MM[16*chan + 14]*Profile->rstdu[i]
			    + RunMode->MM[16*chan + 15]*Profile->rstdv[i]); 

      if(Diagnose){
	if(i==PeakBin && hdr->obs.ChanFreq[chan]==1432.0){
	  dump++;
	  printf("CHANNEL %lf -- DUMP %d\n",hdr->obs.ChanFreq[chan],dump);
	  if(dump == 10 ){
	    if ((fptest = fopen("check_mueller.dat","w")) == NULL) {
	      printf("Cannot open check_mueller.dat. Exiting...\n"); 
	      exit(1);
	    }
	    //	  printf("DUMP %d\n\n",dump++);
	    fprintf(fptest,"S_before =  %9.5f\n",Profile->rstds[i]);
	    fprintf(fptest,"S_after  =  %9.5f + %9.5f + %9.5f + %9.5f = %9.5f\n\n",
		   RunMode->MM[0]*Profile->rstds[i],
		   RunMode->MM[1]*Profile->rstdq[i],
		   RunMode->MM[2]*Profile->rstdu[i],
		   RunMode->MM[3]*Profile->rstdv[i],
		   TempProf.rstds[i]);
	    fprintf(fptest,"Q_before =  %9.5f\n",Profile->rstdq[i]);
	    fprintf(fptest,"Q_after  =  %9.5f + %9.5f + %9.5f + %9.5f = %9.5f\n",
		   RunMode->MM[4]*Profile->rstds[i],
		   RunMode->MM[5]*Profile->rstdq[i],
		   RunMode->MM[6]*Profile->rstdu[i],
		   RunMode->MM[7]*Profile->rstdv[i],
		   TempProf.rstdq[i]);
	    fprintf(fptest,"U_before =  %9.5f\n",Profile->rstdu[i]);
	    fprintf(fptest,"U_after  =  %9.5f + %9.5f + %9.5f + %9.5f = %9.5f\n",
		   RunMode->MM[8]*Profile->rstds[i],
		   RunMode->MM[9]*Profile->rstdq[i],
		   RunMode->MM[10]*Profile->rstdu[i],
		   RunMode->MM[11]*Profile->rstdv[i],
		   TempProf.rstdu[i]);
	    fprintf(fptest,"V_before =  %9.5f\n",Profile->rstdv[i]);
	    fprintf(fptest,"V_after  =  %9.5f + %9.5f + %9.5f + %9.5f = %9.5f\n\n\n",
		   RunMode->MM[12]*Profile->rstds[i],
		   RunMode->MM[13]*Profile->rstdq[i],
		   RunMode->MM[14]*Profile->rstdu[i],
		   RunMode->MM[15]*Profile->rstdv[i],
		   TempProf.rstdv[i]);
	    fflush(stdout);	
	    fclose(fptest);
	  }
	}
      }
      
    }
    
  }
  
  /* Output ThetaBB-corrected file */

  memcpy(Profile,&TempProf,sizeof(struct StdProfs));

  //  if(Diagnose) fclose(fptest);

  return 0;
}




