#include <stdio.h>
#include <math.h>
/* #include "ASPHeader.h" */
#include "fitsio.h"
#include "ASPCommon.h"

int ReadCal(struct ASPHdr *hdr, struct RunVars *RunMode, 
	    struct CalVars *CalMode, double **JyPerCount)
{

  int           i, i_dump, j, chan;
  char          Calfile[128];

  int           NumChansFound, ChanMatch[NCHMAX];

  char          line[128];
  double        CalFreq;

  FILE          *Fcal;


  /* First, grab ThetaBB from file if included in command line */
#if 0
  if (RunMode->ThetaBB){
    if ((Ftheta = fopen(CalMode->Thetafile,"r")) == NULL){
      printf("Could not open phase offset (ThetaBB) calibration file.\n");
      printf("Exiting...\n");fflush(stdout);
      return -10;
    }

    if(RunMode->Verbose){
      printf("\n\nRESULTS FROM CALCULATION OF PHASE SHIFT (THETABB):\n");
      printf(" Freq(MHz)  ThetaBB\n");
      printf(" ---------  -------\n");fflush(stdout);
    }

    chan = 0;
    while ((fgets(line,100,Ftheta) != NULL)){
      
      if(chan == hdr->obs.NChan) break;
      /* Grab ThetaBB values */
      sscanf(line,"%lf%lf",&CalFreq, &ThetaBB[chan]);
      /* Check to make sure that frequencies match */
      if(CalFreq != hdr->obs.ChanFreq[chan]){
	printf("Frequencies do not match between ThetaBB file\n");
	printf("file and header.  Exiting...\n");fflush(stdout);
	return -11;
      }

      if(RunMode->Verbose){
	printf("%10.1lf  %7.2lf\n",
	       hdr->obs.ChanFreq[chan],ThetaBB[chan]);fflush(stdout);
      }

      chan++;

    }

  }
#endif


  if(!RunMode->Cal){  /* if we are not calibrating */
    for (i=0;i<hdr->obs.NChan;i++)
      JyPerCount[i][0] = JyPerCount[i][1] = JyPerCount[i][2] = JyPerCount[i][3]
	= 1.0;
    printf("No cal file given.  Will set cal factors to 1.0.\n");fflush(stdout);
    return 1;
  }




  /* If we've made it this far then we're doing cals, so read in JyPerCount's
   * from file*/

  strcpy(Calfile,  CalMode->Calfile);
  printf("Using calibration factors given in %s.\n\n",Calfile);fflush(stdout);




  if ((Fcal = fopen(CalMode->Calfile,"r")) == NULL){
    printf("Could not open calibration file %s.  Exiting...\n",
	   CalMode->Calfile); fflush(stdout);
    return -6;
  }

  if(RunMode->Verbose){
    printf("\n\nCAL FACTORS (Jy/count):\n");
    printf("=======================\n\n");
    printf(" Freq (MHz)     A^2     B^2     Re(A*B)     Im(A*B) \n");
    printf(" ----------     ---     ---     -------     -------\n\n");
  }

  /* Read file and count to be sure all channels in main data file 
     are accounted for in cal file */
  memset(&ChanMatch[0],0,NCHMAX*sizeof(int));
  NumChansFound = 0;
  chan = 0;
  while ((fgets(line,100,Fcal) != NULL)){
    
    //if(chan == hdr->obs.NChan) break;
    /* Grab cal values */
    sscanf(line,"%lf%lf%lf%lf%lf",
	   &CalFreq, 
	   &JyPerCount[chan][0], &JyPerCount[chan][1],
	   &JyPerCount[chan][2], &JyPerCount[chan][3]);

    for (j=0;j<hdr->obs.NChan;j++){
      //      if(hdr->obs.ChanFreq[j] == CalFreq) {
      /* Test for equality -- do this way since these are doubles */
      if(fabs(hdr->obs.ChanFreq[j] - CalFreq) < DBLEPS) {
	NumChansFound++;
	ChanMatch[j]++;
	CalMode->CalIndex[j] = chan;
	/*	if(RunMode->Verbose){
		printf("%7.1lf   %11.7lf %11.7lf %11.7lf %11.7lf\n",
		hdr->obs.ChanFreq[chan],
		JyPerCount[chan][0],JyPerCount[chan][1],JyPerCount[chan][2],
		JyPerCount[chan][3]);fflush(stdout);	  
		} */
	break;
      }
    }


    chan++;
    
  }

  /* If there is a mismatch between any of the input data ad cal channel sets, 
     then add the extra data channel to the omit matrix by zapping all dumps 
     from a given channel */
  if(NumChansFound != hdr->obs.NChan){
    for (j=0;j<hdr->obs.NChan;j++) {
      if(ChanMatch[j] == 0) {

	for (i_dump=0; i_dump<RunMode->NDumps; i_dump++)
	  RunMode->OmitFlag[i_dump*hdr->obs.NChan + j]=1;
	
	//	if(RunMode->AddChans){
	  RunMode->ZapChan[j]=1;
	  if(j==0){
	    printf("Warning: All channels in input data file are NOT \n");
	    printf("  accounted for in cal file. The following channels\n");
	    printf("   will be omitted:\n\n");
	  }
	  //}
	/*	else{
	  JyPerCount[j][0] = JyPerCount[j][1] 
	    = JyPerCount[j][2] = JyPerCount[j][3] = 1.0;	  
	  printf("Warning: All channels in input data file are NOT \n");
	  printf("  accounted for in cal file. The following channels\n");
	  printf("  will have calibration factor set to 1:\n");	  
	  } */
	printf("Channel %d:  %lf MHz\n",j,hdr->obs.ChanFreq[j]);
	fflush(stdout);
      }
    }
    //    return -4;

    printf("A total of %d channels were omitted from the final calibrated data scan.\n", 
	   hdr->obs.NChan - NumChansFound);
  }

  
  if(RunMode->Verbose) {
    printf("CHANNEL             CAL FACTORS\n");
    for (j=0;j<hdr->obs.NChan;j++) {
      printf("%lf --> %lf   %lf   %lf   %lf\n",hdr->obs.ChanFreq[j],JyPerCount[j][0],JyPerCount[j][1],JyPerCount[j][2], JyPerCount[j][3]);
    }
  }    



  return 0;

}
