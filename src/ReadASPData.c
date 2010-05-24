#include <stdio.h>
#include <math.h>
/* #include "fitsio.h" */
#include "ASPCommon.h"

int ReadASPData(struct ASPHdr *hdr, struct SubHdr *subhdr, 
		 struct RunVars *RunMode, fitsfile *Fin, int nscan, 
		 int NPtsProf, double **ASquared, double **BSquared, 
		 double **ReAconjB, double **ImAconjB, int **SampleCount,
		 char **HeadLine)
{

  int  i, j, status=0, datacolnum, refcolnum, anynull;
  int  NColumns;
  
  char Rawfile[100];

  FILE *s;


  fits_get_num_cols(Fin, &NColumns, &status);status=0;


  /* malloc array for each channel based on profile size */
  for(i=0;i<hdr->obs.NChan;i++){
    ASquared[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ASquared[i][0], NPtsProf);
    BSquared[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&BSquared[i][0], NPtsProf);
    ReAconjB[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ReAconjB[i][0], NPtsProf);
    ImAconjB[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ImAconjB[i][0], NPtsProf);
    SampleCount[i] = (int   *)malloc(NPtsProf*sizeof(int));
    IZero(&SampleCount[i][0], NPtsProf);
  }

  /* get middle dump time  */
  if(!strcmp(hdr->gen.HdrVer,"Ver1.0")) 
    fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", &(subhdr->DumpMiddleSecs), 
		  NULL, &status); status = 0;
  
  /* move to ref table for newer version */
  if(!strcmp(hdr->gen.HdrVer,"Ver1.0.1")) {
    fits_read_key(Fin, TDOUBLE, "MIDSECS", &(subhdr->DumpMiddleSecs), 
		  NULL, &status); status = 0;
    refcolnum = 0;
    fits_read_col(Fin, TDOUBLE, ++refcolnum, 1, 1, (long)hdr->obs.NChan, NULL, 
		  subhdr->DumpRefPhase,  &anynull, &status); 
    fits_read_col(Fin, TDOUBLE, ++refcolnum, 1, 1, (long)hdr->obs.NChan, NULL, 
		  subhdr->DumpRefPeriod, &anynull, &status); 

   if(RunMode->Verbose){
      printf("Dump %d:  TIME OF DUMP = %lf\n",nscan,subhdr->DumpMiddleSecs); 
      printf("          CHANNEL (MHz)   REF. PHASE   REF. PERIOD (s)\n");
      printf("          -------------   ----------   ---------------\n");
      for(i=0;i<hdr->obs.NChan;i++) 
	printf("          %13.1lf%13.8lf%18.11lf\n", hdr->obs.ChanFreq[i],
	       subhdr->DumpRefPhase[i], subhdr->DumpRefPeriod[i]);       
      fflush(stdout);
    }

    /* Move to next HDU, i.e. data table for this dump */
    fits_movrel_hdu(Fin, 1, NULL, &status);
  }

  datacolnum = 0;

  for(i=0;i<hdr->obs.NChan;i++) {

    fits_read_col(Fin, TDOUBLE,  ++datacolnum, 1, 1, NPtsProf, NULL, 
		  &ASquared[i][0],    &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++datacolnum, 1, 1, NPtsProf, NULL, 
		  &BSquared[i][0],    &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++datacolnum, 1, 1, NPtsProf, NULL, 
		  &ReAconjB[i][0],    &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++datacolnum, 1, 1, NPtsProf, NULL, 
		  &ImAconjB[i][0],    &anynull, &status); 
    if(!RunMode->OldFits)
      fits_read_col(Fin, TLONG,  ++datacolnum, 1, 1, NPtsProf, NULL, 
		    &SampleCount[i][0], &anynull, &status); 


/* Normalize */
    
    if(!RunMode->OldFits){
      for(j=0;j<NPtsProf;j++){
        ASquared[i][j] /= (double)SampleCount[i][j];
        BSquared[i][j] /= (double)SampleCount[i][j];
        ReAconjB[i][j] /= (double)SampleCount[i][j];
        ImAconjB[i][j] /= (double)SampleCount[i][j];

      }
    }

  }


  for(i=0;i<hdr->obs.NChan;i++) {

   /* read in reference phase, period if old version */


    if(!strcmp(hdr->gen.HdrVer,"Ver1.0")){
      fits_read_key(Fin, TDOUBLE, "DUMPREFPER",   &(subhdr->DumpRefPeriod[i]),  
		    NULL, &status); status = 0;
      fits_read_key(Fin, TDOUBLE, "DUMPREFPHASE", &(subhdr->DumpRefPhase[i]),   
		    NULL, &status); status = 0; 
                                                                              
      if(RunMode->Verbose && i==0){
	printf("Dump %d:  TIME OF DUMP = %lf\n",nscan, subhdr->DumpMiddleSecs);
	printf("          REF. PHASE = %lf, REF. PERIOD = %lf s\n", 
	       subhdr->DumpRefPhase[i], subhdr->DumpRefPeriod[i]);fflush(stdout);
      }
    }



   sprintf(HeadLine[i],"# %.1f %.7f %.10f %ld %.3f %.3f %d %s %d %s %.10f",
	    (double)hdr->obs.IMJDStart, subhdr->DumpMiddleSecs, 
	    subhdr->DumpRefPeriod[i], (long)1,hdr->obs.ChanFreq[i], hdr->obs.DM, 
	    RunMode->NBinsOut, hdr->obs.ObsvtyCode, 1, hdr->target.PSRName, 
	    subhdr->DumpRefPhase[i]);

    /* Open and write raw  file if asked for by user */
    if (RunMode->MakeRaw){
      sprintf(Rawfile,"%s.%4.4d.%4.4d.prof.raw.asc",RunMode->OutfileRoot,
	      (int)(hdr->obs.ChanFreq[i]),nscan);
      s = fopen(Rawfile, "w");
      if (s==NULL) {
	printf("ASPFitsReader: Can't open file %s\n",Rawfile);
	fflush(stdout);
	exit(1);
      }

      fprintf(s,"%s\n",HeadLine[i]);

      for (j=0; j<NPtsProf; j++) {
 	fprintf(s,"%d  %lf  %lf  %lf  %lf  %d\n",j, 
		ASquared[i][j], BSquared[i][j], ReAconjB[i][j], ImAconjB[i][j],
		SampleCount[i][j]); 
      }
      fclose(s);

    }

  }

  return 1;

}
