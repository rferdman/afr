#include <stdio.h>
#include <math.h>
/* #include "fitsio.h" */
#include "ASPCommon.h"
#include "phase_calc.h"

int ReadPSRFITSData(struct ASPHdr *hdr, struct SubHdr *subhdr, 
		    struct RunVars *RunMode, fitsfile *Fin, 
		    int i_dump, int NPtsProf, 
		    double **ASquared, double **BSquared, 
		    double **ReAconjB, double **ImAconjB)
{
  
  int  i_chan, i_poly, i_poln, i_bin, i_data;
  int  hdutype, anynul, status=0;
  int  colnum, colnum_wts, colnum_scl, colnum_offs, colnum_data;
  int  NColumns, n_poly, FileStartSecs;
  static int  n_suboffs=0;   
  double  *tempprof;
  //short int  *tempprof;
  float  *Weights, *Offsets, *Scales;
  static double TotalTDump=0.;
  double ref_mjd, DumpMiddleDays, RefPhase, RefFreq;
  char Rawfile[100], tempstr[8];
  struct Polyco Polycos[MAX_PC_SETS];

  FILE *s;


  fits_get_num_cols(Fin, &NColumns, &status);status=0;


  /* malloc array for each channel based on profile size */
  for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
    if (i_dump==0) {
      ASquared[i_chan]=(double *)malloc(hdr->redn.RNBinTimeDump*sizeof(double));
      BSquared[i_chan]=(double *)malloc(hdr->redn.RNBinTimeDump*sizeof(double));
      ReAconjB[i_chan]=(double *)malloc(hdr->redn.RNBinTimeDump*sizeof(double));
      ImAconjB[i_chan]=(double *)malloc(hdr->redn.RNBinTimeDump*sizeof(double));
    }
    /*    SampleCount[i_chan] = (int    *)malloc(NPtsProf*sizeof(int));
	  IZero(&SampleCount[i_chan][0], NPtsProf); */
    DZero(&ASquared[i_chan][0], hdr->redn.RNBinTimeDump);
    DZero(&BSquared[i_chan][0], hdr->redn.RNBinTimeDump);
    DZero(&ReAconjB[i_chan][0], hdr->redn.RNBinTimeDump);
    DZero(&ImAconjB[i_chan][0], hdr->redn.RNBinTimeDump);
  }
  

  /* malloc temporary arrays for profile data */
  ////////tempprof = (int *)malloc(NPtsProf*sizeof(int));
  ////////tempprof = (double *)malloc(hdr->redn.RNBinTimeDump*hdr->obs.NChan*hdr->obs.NPoln*sizeof(double));
  /*   tempprof = (double *)malloc(hdr->redn.RNBinTimeDump*sizeof(double)); */
  tempprof = (double *)malloc(hdr->redn.RNBinTimeDump*hdr->obs.NPoln*sizeof(double));
  ////////tempprof = (short int *)malloc(hdr->redn.RNBinTimeDump*sizeof(short int));
  /*  tempA = (int *)malloc(NPtsProf*sizeof(int));
  tempB = (int *)malloc(NPtsProf*sizeof(int));
  tempReAconjB = (int *)malloc(NPtsProf*sizeof(int));
  tempImAconjB = (int *)malloc(NPtsProf*sizeof(int)); */

  /* malloc scaling, offset, and weighting values */
  Weights = (float *)malloc(hdr->obs.NChan*sizeof(float));
  Offsets = (float *)malloc(hdr->obs.NPoln*hdr->obs.NChan*sizeof(float));
  Scales  = (float *)malloc(hdr->obs.NPoln*hdr->obs.NChan*sizeof(float));
  if(Scales == NULL || Offsets == NULL || Weights == NULL) {
    fprintf(stderr, "ReadPSRFITSData ERROR: Memory allocation failed\n");
    return -1;
  }

  /* Read in polycos from FITS file -- first move to POLYCO extension */
  if(fits_movnam_hdu(Fin, BINARY_TBL, "POLYCO", 0, &status)) {
    fprintf(stderr, "ReadPSRFITSData WARNING: POLYCO table does not exist\n");
    fprintf(stderr, "Continuing, with no pulse phase information.\n");
    /* Reset status = 0 so that further fitsio routines will work */
    status = 0;
    /* Assign reference frequency to be centre frequency for the purposes of 
     initial alignment */
    n_poly = 0;
    Polycos[0].FSkyRef = hdr->obs.FSkyCent;
  }
  /* Read in number of rows, ie. number of polyco sets */
  else{
    if(fits_read_key(Fin, TINT, "NAXIS2", &n_poly, NULL, &status)) {
      fprintf(stderr, "ReadPSRFITSData ERROR:  Could not get NAXIS2 from POLYCO extension.");
      return -1;
    }
    /* Now read in actual information for each polyco set in data file */
    /* NSPAN -- Valid span of polyco set in minutes */
    for(i_poly=0; i_poly<n_poly; i_poly++){
      if(fits_get_colnum (Fin, CASEINSEN, "NSPAN", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No NSPAN in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TINT, colnum, 1+i_poly, 1, 1, NULL, &Polycos[i_poly].NMinutes, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read NSPAN...\n");
	return -1;
      }
      
      /* NCOEFF -- Number of coefficients*/
      if(fits_get_colnum (Fin, CASEINSEN, "NCOEF", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No NCOEF in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TINT, colnum, 1+i_poly, 1, 1, NULL, &Polycos[i_poly].NCoeff, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read NCOEF...\n");
	return -1;
      }
      
      /* REF_MJD -- Reference MJD for this polyco set */
      if(fits_get_colnum (Fin, CASEINSEN, "REF_MJD", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No REF_MJD in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, &ref_mjd, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read REF_MJD...\n");
	return -1;
      }
      Polycos[i_poly].MjdMidInt = floor(ref_mjd);
      Polycos[i_poly].MjdMidFrac = ref_mjd - floor(ref_mjd);
      
      printf("ref_mjd:                %lf\n", ref_mjd);
      printf("Polycos:  MjdMidInt:    %lf\n",  Polycos[i_poly].MjdMidInt);
      printf("Polycos:  MjdMidFrac:   %lf\n\n",  Polycos[i_poly].MjdMidFrac);
      
      /* REF_FREQ -- Reference observing frequency for this polyco set */
      if(fits_get_colnum (Fin, CASEINSEN, "REF_FREQ", &colnum, &status)){
	fprintf(stderr,"ReadPSRFITSData ERROR: No REF_FREQ in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, &Polycos[i_poly].FSkyRef, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read REF_FREQ...\n");
	return -1;
      }
      
      /* REF_PHS -- Reference phase for this polyco set */
      if(fits_get_colnum (Fin, CASEINSEN, "REF_PHS", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No REF_PHS in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, &Polycos[i_poly].PhRotRef, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read REF_PHS...\n");
	return -1;
      }
      
      /* REF_F0 -- Reference rotation frequency for this polyco set */
      if(fits_get_colnum (Fin, CASEINSEN, "REF_F0", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No REF_F0 in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, &Polycos[i_poly].FRotRef, &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read REF_F0...\n");
	return -1;
      }
      
      /* COEFF -- Actual polynomial coefficients for this set */
      if(fits_get_colnum (Fin, CASEINSEN, "COEFF", &colnum, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: No COEFF in POLYCO table.\n");
	return -1;
      }
      if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, Polycos[i_poly].NCoeff, NULL, &Polycos[i_poly].Coeff[0], &anynul, &status)){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read REF_F0...\n");
	return -1;
      }
    }
  }
  
  
  
  /* Move to SUBINT extension */
  if(fits_movnam_hdu(Fin, BINARY_TBL, "SUBINT", 0, &status)) {
    fprintf(stderr, "ReadPSRFITSData ERROR: SUBINT extension does not exist!\n");
    return -1;
  } 

  /* Read in number of rows before this file if it is a file in a series, so 
     that we can add the offset to the number of subints in this file */
  if (i_dump==0) { // first scan only
    if(fits_read_key(Fin, TINT, "NSUBOFFS", &n_suboffs, NULL, &status)) {
      if(!strncmp(hdr->gen.ObsMode, "CAL", 3)){
	fprintf(stderr, "ReadPSRFITSData WARNING:  Could not get NSUBOFFS from SUBINT extension.\n");
	printf("Proceeding with CAL file reading.\n");	
	status=0;  /* Otherwise next fits routine will crash */

      } 
      else {
	fprintf(stderr, "ReadPSRFITSData WARNING:  Could not read NSUBOFFS keyword from SUBINT extension.\n");
	printf("Proceeding with reading of data.  Time stamps likely not correct.\n");	
	n_suboffs=0;
	//sleep(2);
	status=0;  /* Otherwise next fits routine will crash */
	//  fits_report_error(stderr, status); /* print any error message */
	
      }
    }

  }
  
/* Length of integration */
  if(fits_get_colnum (Fin, CASEINSEN, "TSUBINT", &colnum, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No TSUBINT in FITS file?\n");
    return -1;
  }
  if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_dump, 1, 1, NULL, &hdr->redn.TDump, &anynul, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read TSUBINT...\n");
    return -1;
  }
  /* If first subint in file, start by adding n_suboffs' worth of scan time 
     to total TDump -- needed to get correct time of each subint in current 
     data file -- and assume TDump for all dumps is the same as this first 
     dump in the current file */
  if(i_dump==0) {
    TotalTDump = ((double)n_suboffs)*hdr->redn.TDump;
    /* Get the number of seconds after UT 0:00 this particular scan began */
    FileStartSecs = hdr->obs.StartTime + 
      (int)floor((double)n_suboffs*hdr->redn.TDump);
    /* Now with all the required info, we can construct the 
       output file name */
    sprintf(RunMode->OutfileRoot, "%s.%d.%d.%s.%s", 
	    hdr->target.PSRName, 
	    hdr->obs.IMJDStart + (int)(FileStartSecs/86400), 
	    FileStartSecs%86400, 
	    hdr->obs.ObsvtyCode, hdr->gen.BEName);
  
  }
  TotalTDump += hdr->redn.TDump;


  /* Now calculate the mid-scan dump time based on TotalTDump, TDump, Start time, and current dump number */
  
  /* Check that polyco falls inside NMinutes/2 of subint time */
  
  subhdr->DumpMiddleSecs =  (double)hdr->obs.StartTime + 
    TotalTDump - hdr->redn.TDump/2.;

  /* Need to recalculate RefPhase and RefPeriod based on polycos or something */
  DumpMiddleDays = subhdr->DumpMiddleSecs /86400.;
  
  printf("\n\nDUMP %d:\n", i_dump);
  printf("========\n");
  printf("Dump length:                             %lf seconds\n", 
	 hdr->redn.TDump);
  printf("Time since scan start:                   %lf seconds\n", 
	 TotalTDump);
  printf("Midpoint dump time (secs since start):   %lf\n\n", 
	 subhdr->DumpMiddleSecs);
  printf("Midpoint dump time (days since start):   %lf\n", 
	 DumpMiddleDays);
  printf("Dump time (MJD):   %lf\n", 
	 hdr->obs.IMJDStart+DumpMiddleDays);


  if(PhaseCalc(Polycos, n_poly, hdr->obs.IMJDStart, DumpMiddleDays, 
	       &RefPhase, &RefFreq) < 0){
    printf("ReadPSRFITSData WARNING: Could not find any polyco sets to \n");
    printf("match MJD for PSRFITS file.\n");
    printf("Phase and period information not known.  Will not be able to \n");
    printf("realign using polycos.\n\n");
    //return -1;
  }

  /*** IF NO VALID POLYCOS BECAUSE FITS FILE DATES ARE OFF... 
       THEN FIND PHASE WITH FFTFIT AND SET REFPERIOD TO 0 ***/
  /*** AND IF NO NEW POLYCO FILE IS GIVEN ON COMMAND LINE, 
       GIVE A WARNING THAT IT'S PROBABLY NEEDED ***/


  /* In general, PSRFITS data does not store different polycos for different 
     observing frequency channels, so will have to assign each period/phase 
     combo by hand here. Can correct by inputting a new polyco file which will 
     be used later in data processing chain */
  for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
    subhdr->DumpRefPeriod[i_chan] = 1./RefFreq;
    subhdr->DumpRefPhase[i_chan]  = RefPhase;
  }

 

  /* Find colnum for each of weights, offsets, and scales */
  if(fits_get_colnum (Fin, CASEINSEN, "DAT_WTS", &colnum_wts, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No data weights in fits file?\n");
    return -1;
  }
  if(fits_read_col(Fin, TFLOAT, colnum_wts, 1+i_dump, 1, hdr->obs.NChan, NULL, &Weights[0], &anynul, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read weights...\n");
    return -1;
  }

  
  if(fits_get_colnum (Fin, CASEINSEN, "DAT_OFFS", &colnum_offs, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No data offsets in fits file?\n");
    return -1;
  }
  
  if(fits_get_colnum (Fin, CASEINSEN, "DAT_SCL", &colnum_scl, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No data scaling in fits file?\n");
    return -1;
  }
  
  if(fits_get_colnum (Fin, CASEINSEN, "DATA", &colnum_data, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No data in fits file?\n");
    return -1;
  }

  for(i_poln=0; i_poln<hdr->obs.NPoln; i_poln++){

    //    if(fits_read_col(Fin, TFLOAT, colnum_offs, 1+i_dump, 1, hdr->obs.NPoln*hdr->obs.NChan, NULL, &Offsets[0], &anynul, &status)){
    if(fits_read_col(Fin, TFLOAT, colnum_offs, 1+i_dump, 1+i_poln, hdr->obs.NChan, NULL, &Offsets[i_poln*hdr->obs.NChan], &anynul, &status)){
      fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read offsets...\n");
      return -1;
    }
    
    //    if(fits_read_col(Fin, TFLOAT, colnum_scl, 1+i_dump, 1, hdr->obs.NPoln*hdr->obs.NChan, NULL, &Scales[0], &anynul, &status)){
    if(fits_read_col(Fin, TFLOAT, colnum_scl, 1+i_dump, 1+i_poln, hdr->obs.NChan, NULL, &Scales[i_poln*hdr->obs.NChan], &anynul, &status)){
      fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read scales...\n");
      return -1;
    }
    
  }

#if 0
  /*  Check Offsets and Scales for polarization 0 */
  /* for (i_chan=0; i_chan< hdr->obs.NChan; i_chan++)
    printf ("%.2f   ", Offsets[i_chan]);
    printf("\n\n"); */
  for (i_chan=0; i_chan< hdr->obs.NChan; i_chan++)
    printf ("%.2f   ", Scales[i_chan]);
  printf("\n\n");
  sleep(10);
#endif


  int initflag = 0;
  int typecode = 0;
  long repeat = 0;
  long width = 0;
  status=0;
  i_data = 1;

  //printf("\nnbin = %d\npoln = %d\nnchan = %d\n\n", hdr->redn.RNBinTimeDump,hdr->obs.NPoln, hdr->obs.NChan );

  fits_get_coltype (Fin, colnum_data, &typecode, &repeat, &width, &status);  
  if(typecode == TSHORT){
    //free(tempprof);
    //    short int *tempprof;
    //tempprof = (short int *)malloc(hdr->redn.RNBinTimeDump*sizeof(short int));
    printf("DATA type is TSHORT\n");
  }
  else if(typecode == TFLOAT){
    // free(tempprof);
    // float *tempprof;
    // tempprof = (float *)malloc(hdr->redn.RNBinTimeDump*sizeof(float));
    printf("DATA type is TFLOAT\n");    
  }
  else if(typecode == TINT){
    //free(tempprof);
    //int *tempprof;
    //tempprof = (int *)malloc(hdr->redn.RNBinTimeDump*sizeof(int));
    printf("DATA type is TINT\n");
  }
  else{
    printf("DATA type %d is unknown.\n",typecode);
  }
  // exit(0);

      //i_poln=0; // start at ASquared
  
  
      
  for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++) {
    for(i_poln=0; i_poln<hdr->obs.NPoln; i_poln++) {
      
      i_data = 1 + i_poln*hdr->redn.RNBinTimeDump*hdr->obs.NChan + 
	i_chan*hdr->redn.RNBinTimeDump;
      
      fits_read_col(Fin, TDOUBLE, colnum_data, 1+i_dump, i_data, 
		    hdr->redn.RNBinTimeDump, 
		    NULL, &(tempprof[i_poln*hdr->redn.RNBinTimeDump]), 
		    &initflag, &status);
      //    if(fits_read_col(Fin, TINT, colnum_data, 1+i_dump, 1+i_poln*hdr->redn.RNBinTimeDump*hdr->obs.NChan + i_chan*hdr->redn.RNBinTimeDump, hdr->redn.RNBinTimeDump, NULL, &tempprof[0], &anynul, &status)){
      //      i_data += hdr->redn.RNBinTimeDump;            
      
      if(status != 0){
	fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read profile...\n");
	return -1;
      }
      
    }
 

    /* Convert Stokes parameters back to cross-products if the PSRFITS data 
       is in IQUV format already.  Otherwise, leave alone: */
    /* Check that POL_TYPE is "IQUV" or "STOKE", etc. */
    
    if (!strncasecmp(&hdr->obs.PolnType[0], "IQUV", 4) ||
	!strncasecmp(&hdr->obs.PolnType[0], "STOK", 4)) {
      if(!strcmp("L",hdr->gen.FEPol)){  /* dual linear */
	if(RunMode->Swap){
	  for(i_bin=0; i_bin<hdr->redn.RNBinTimeDump; i_bin++) {
	    /* If LIN, with A/B pol swap */
	    ASquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] -
	       Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]) ; 
	    
	    BSquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] +
	       Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]) ; 
	    
	    ReAconjB[i_chan][i_bin] = (double)
	      (Scales[2 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[2 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[2 * hdr->obs.NChan + i_chan]);
	    
	    ImAconjB[i_chan][i_bin] = -1.0 * (double)
	      (Scales[3 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[3 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[3 * hdr->obs.NChan + i_chan]);
	  }
	}
	
	else {
	  for(i_bin=0; i_bin<hdr->redn.RNBinTimeDump; i_bin++) {
	    /* If LIN, no swap */
	    ASquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] +
	       Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]) ; 
	    
	    BSquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] -
	       Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]) ; 
	    
	    ReAconjB[i_chan][i_bin] = (double)
	      (Scales[2 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[2 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[2 * hdr->obs.NChan + i_chan]);
	    
	    ImAconjB[i_chan][i_bin] = (double)
	      (Scales[3 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[3 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[3 * hdr->obs.NChan + i_chan]);
	  }
	}
      }
      else if (!strcmp("C",hdr->gen.FEPol)) { /* dual circular */
	
	if(RunMode->Swap){
	  for(i_bin=0; i_bin<hdr->redn.RNBinTimeDump; i_bin++) {
	    /* If CIRC, with A/B pol swap */
	    ASquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] -
	       Scales[4 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[4 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[4 * hdr->obs.NChan + i_chan]) ; 
	    
	    BSquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] +
	       Scales[4 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[4 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[4 * hdr->obs.NChan + i_chan]) ; 
	    
	    ReAconjB[i_chan][i_bin] = (double)
	      (Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]);
	    
	    ImAconjB[i_chan][i_bin] = -1.0 * (double)
	      (Scales[2 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[2 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[2 * hdr->obs.NChan + i_chan]);
	  }
	}
	else {
	  for(i_bin=0; i_bin<hdr->redn.RNBinTimeDump; i_bin++) {
	    /* If CIRC, no swap */
	    ASquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] +
	       Scales[4 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[4 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[4 * hdr->obs.NChan + i_chan]) ; 
	    
	    BSquared[i_chan][i_bin] = (double) 
	      (Scales[0 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[0 * hdr->obs.NChan + i_chan] -
	       Scales[4 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[4 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[4 * hdr->obs.NChan + i_chan]) ; 
	    
	    ReAconjB[i_chan][i_bin] = (double)
	      (Scales[1 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[1 * hdr->obs.NChan + i_chan]);
	    
	    ImAconjB[i_chan][i_bin] = (double)
	      (Scales[2 * hdr->obs.NChan + i_chan]*
	       ((float)tempprof[2 * hdr->redn.RNBinTimeDump + i_bin]) +
	       Offsets[2 * hdr->obs.NChan + i_chan]);
	  }
	}
	
      }
      else {  /* Neither Linear nor Circular polarizations */
	fprintf(stderr, "ReadPSRFITSData ERROR: Unrecognized polarization ");
	fprintf(stderr, "basis -- %s\n.", hdr->gen.FEPol);
	return -1;
      }
    }
    /* If not in Stokes format, then for now assume AA-BB-Re-Im format */
    else{
      
      for(i_bin=0; i_bin<hdr->redn.RNBinTimeDump; i_bin++){
	
	// i_data = i_chan*hdr->redn.RNBinTimeDump*hdr->obs.NPoln + i_poln*hdr->redn.RNBinTimeDump + i_bin;            
	//	i_data = i_poln*hdr->redn.RNBinTimeDump*hdr->obs.NChan + i_chan*hdr->redn.RNBinTimeDump + i_bin;       
	
	//   ASquared[i_chan][i_bin] = (double)(Weights[i_chan]*
	//	if(i_poln==0) 
	
	/* Now factor in weighting, scaling, and offsets:
	   weight*(data*scale + offsets) */
	
	
	ASquared[i_chan][i_bin] =  (double) //tempprof[i_data+i_bin];
	  (Scales[0 * hdr->obs.NChan + i_chan]*
	   ((float)tempprof[0 * hdr->redn.RNBinTimeDump + i_bin]) +
	   Offsets[0 * hdr->obs.NChan + i_chan]);  	
	//	if(i_poln==1)
	BSquared[i_chan][i_bin] = (double)
	  (Scales[1 * hdr->obs.NChan + i_chan]*
	   ((float)tempprof[1 * hdr->redn.RNBinTimeDump + i_bin]) +
	   Offsets[1 * hdr->obs.NChan + i_chan]);
	//	if(i_poln==2)
	ReAconjB[i_chan][i_bin] = (double)
	  (Scales[2 * hdr->obs.NChan + i_chan]*
	   ((float)tempprof[2 * hdr->redn.RNBinTimeDump + i_bin]) +
	   Offsets[2 * hdr->obs.NChan + i_chan]);
	//	if(i_poln==3)
	ImAconjB[i_chan][i_bin] = (double)
	  (Scales[3 * hdr->obs.NChan + i_chan]*
	   ((float)tempprof[3 * hdr->redn.RNBinTimeDump + i_bin]) +
	   Offsets[3 * hdr->obs.NChan + i_chan]);
      }
    }
    //printf("\n");
    
   

    
  }
  
  
// printf("\n\nnbin = %d\npoln = %d\nnchan = %d\n", hdr->redn.RNBinTimeDump,hdr->obs.NPoln, hdr->obs.NChan );

      //exit(0);
      
      

      
      /*   
	   sprintf(HeadLine[i],"# %.1f %.7f %.10f %ld %.3f %.3f %d %s %d %s %.10f",
	   (double)hdr->obs.IMJDStart, subhdr->DumpMiddleSecs, 
	   subhdr->DumpRefPeriod[i], (long)1,hdr->obs.ChanFreq[i], hdr->obs.DM, 
	   RunMode->NBinsOut, hdr->obs.ObsvtyCode, 1, hdr->target.PSRName, 
       subhdr->DumpRefPhase[i]);
       
  */


  /* Free up some memory */
  /*  free(tempA);
  free(tempB);
  free(tempReAconjB);
  free(tempImAconjB); */
  free(tempprof);
  free(Weights);
  free(Offsets);
  free(Scales);
  
  
  
  
  return 0;
}
