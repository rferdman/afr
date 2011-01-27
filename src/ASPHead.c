/* Utility to print out primary header information from fits file */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPHeader.h"
#include "fitsio.h"
#include "HeadCmdLine.h"

int main(int argc, char **argv)
{

  int      i, j, indx, n_file, special_key, knul=0, anynull, hdutype, status=0;
  int      LastSlashIndex;
  int      NDump, NChan, NumHDU;
  long     nrows=0;
  double   ChanFreq, BW, CentFreq, DM;
  double   TDump, StartTime, DumpMiddleSecs, ScanLen;
  char     KeyValue[32], err_text[256];
  char     cur_file[256], cur_file_short[64];
  char     HdrVer[16];
  char     aspoutstr[11], dumprefstr[11];
  fitsfile *Fin;
  Cmdline  *Cmd;


  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /*  InfileP = (argc > 1); printf("InfileP = %d\n", InfileP);
  InfileC = (argc - 1); printf("InfileC = %d\n", InfileC);
  Infile = (char **)malloc(InfileC);
  for (i=0; i<InfileC; i++) {
    Infile[i] = (char *)malloc(256);
    if (strlen(argv[i+1]) > 255) {
      fprintf(stderr, "Error: argument %d is too long.\n", i+1);
      exit(1);
    }
    else{
      strcpy(Infile[i], argv[i+1]);
    }
    printf("Argument %d:  %s = %s\n", i+1, Infile[i], argv[i+1]);
    } */
  
  /************* Keyword list output  **************/

  if (Cmd->ListP) {
    if (argc > 2) {
      fprintf(stderr,"\n-list option must be used alone.  Exiting...\n");
      exit(1);
    }
    else{
      printf("\nSome available keywords/info options useable with the ");
      printf("-key flag:\n");
      printf("\nObserving frequency-related options:\n");
      printf("\nNCHAN:       No. of channels\n");
      printf("CFRQ:        Centre frequecy\n");
      printf("BW:          Bandwidth\n");
      printf("STARTCHAN:   First (0th) frequency channel\n");
      printf("ENDCHAN:     Last (NCHAN-th) frequency channel\n");
      printf("ALLCHAN:     List all channels (ALLCHAN keyword output does\n");
      printf("             not follow convention of -col flag)\n");
     
      printf("\nIntegration-related options:\n");
      printf("\nNDUMP:       No. of integrations\n");
      printf("TDUMP:       Individual integration time (s)\n");
      printf("SCANLEN:     Full scan length of input file (s)\n");
      printf("NPTSPROF:    No. of bins in integrated profiles\n");

      printf("\nOther options:\n");

      printf("\nSRC_NAME:    Source Name\n");
      printf("DM:          Dispersion measure\n");
      printf("OBS_MODE:    PSR/CAL -- normal or calibration scan\n");
      printf("OBSVTY:      Observatory code\n");
      printf("FD_POLN:     Polarisation basis (e.g. L or C for linear or\n");
      printf("             circular, respectively)\n");
      printf("RA:          Right ascension of source\n");
      printf("DEC:         Declination of source\n");
      
      printf("\nNote that all -key arguments are NOT case-sensitive.\n");
    }
  }
  else if(Cmd->TDiffP) {
    if (Cmd->KeywordP) {
      fprintf(stderr,"\n-tdiff option must be used alone.  Exiting...\n");
      exit(1);
    }
  }
  else{
    if(!Cmd->KeywordP){
      printf("Must use -key option with appropriate arguments if not\n");
      printf("using -list or -tdiff flag\n");
      exit(1);
    }
    if(!Cmd->InfileP){
      printf("Must use -infile option with input file arguments if not\n");
      printf("using -list or -tdiff flag\n");
      exit(1);
    }
  }


  /*************************************************/


  /* Print header line */
  if(Cmd->KeywordP && Cmd->ColumnP) {
    printf("%s  ","File");
    for(i=0;i<Cmd->KeywordC;i++)
      printf("%s  ",Cmd->Keyword[i]);
    printf("\n");
  }
  else if(Cmd->TDiffP){  // for Tdiff
    printf("File                        Source     MJD              TDIFF [ns]\n");
  }



  for(n_file=0;n_file<Cmd->InfileC;n_file++){
    strcpy(cur_file,Cmd->Infile[n_file]);
    /* Open fits file */
    if(fits_open_file(&Fin, cur_file, READONLY, &status)){
      printf("\nError opening FITS file %s\n",
	     cur_file);
      fits_get_errstatus(status, err_text);
      printf("FITS Error status %d: %s\n",status,err_text);
      exit(1);
    }
    status=0;

    /* get rid of everthing before and including last forward slash 
       in the filename for printing to screen: */
    LastSlashIndex = -1;
    for(i=strlen(cur_file)-1;i>0;i--){
      if(!strncasecmp(&cur_file[i],"/",1)){
	LastSlashIndex = i;
	break;
      }
    }

    strncpy(cur_file_short,&cur_file[LastSlashIndex+1],strlen(cur_file)-LastSlashIndex);


    if(Cmd->TDiffP){

      /******************** TDIFF ***********************/
      
      int    i_dump;
      int    IMJDStart, SMJDStart;
      double StartMJD, MJD, TDiff;
      char   ObsCode[4], psr[12];
      
      /* Check this is Nancay data */
      fits_read_key(Fin, TSTRING, "OBSVTY", ObsCode, NULL, &status); 
      if(strncmp(ObsCode,"f",1)) {
	printf("\nError: File %s does not contain Nancay data (obs code = %s).",
	       cur_file_short,ObsCode);
	printf("  Cannot extract TDIFF from this file.  Exiting.\n");
	exit(2);
      }

      /* Get pulsar name */
      fits_read_key(Fin, TSTRING, "SRC_NAME", psr, NULL, &status); 
   

      /* Calculate start MJD */
      fits_read_key(Fin, TINT, "STT_IMJD", 
		    &(IMJDStart),NULL, &status); status = 0;
      fits_read_key(Fin, TINT, "STT_SMJD", 
		    &(SMJDStart),NULL, &status); status = 0;

      StartMJD = (double)IMJDStart + (double)SMJDStart/86400.;

      /* Find nuber of dumps in this file */
      fits_get_num_hdus(Fin, &NumHDU, &status);
      NDump = NumHDU-3;  // Nancay data is always Ver1.0
      
      /* Loop through each dump and pick out time stamp and phase */
      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT0", 0, &status);
    
      for (i_dump=0;i_dump<NDump;i_dump++){
	fits_movrel_hdu(Fin, 1, NULL, &status);
	fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", 
		      &(DumpMiddleSecs),NULL, &status); status = 0;
	fits_read_key(Fin, TDOUBLE, "DUMPTDIFF", 
		      &(TDiff),NULL, &status); status = 0;
	
	MJD = StartMJD + DumpMiddleSecs/86400.;
	
	printf("%s  ",cur_file_short);
	printf("%s  ",psr);
	printf("%14.8lf   ",MJD);
	printf("%7.3lf   \n",TDiff);
      }

    
    }
    // else





    if(Cmd->ColumnP) printf("%s  ",cur_file_short);
    else printf("\nKeywords in %s:\n\n", cur_file_short);
    
    /* For each keyword requested... */
    for(i=0;i<Cmd->KeywordC;i++){
      special_key=0;
     
      /* Clear Keyword */
      strcpy(KeyValue,"\0");

      /**************** Channel-based keywords: *****************/
      if(!strcasecmp(Cmd->Keyword[i],"NChan") ||      // # of channels
	 !strncasecmp(Cmd->Keyword[i],"CFRQ",4) ||    // centre frequency
	 !strncasecmp(Cmd->Keyword[i],"BW",2) ||      // bandwidth in MHz
	 !strcasecmp(Cmd->Keyword[i],"StartChan") ||  // first channel (MHz)
	 !strcasecmp(Cmd->Keyword[i],"EndChan") ||    // last channel (MHz)
	 !strcasecmp(Cmd->Keyword[i],"AllChan")){     // all channels (MHz)
	/* move to channel table */
	fits_movnam_hdu(Fin, ASCII_TBL, "BECONFIG", 0, &status);
	ffgnrw(Fin, &nrows, &status);

	indx = 1;
	fits_read_col(Fin, TINT, indx++, 1, 1, 1, &knul,
		      &NChan, &anynull, &status); 

	/* Write heading for AllChans keyword if requested */
	if (!strcasecmp(Cmd->Keyword[i],"AllChan")) 
	  printf("\nAll Channels:\n");
	/* Read NChans number, copy to keyword */
	if (!strcasecmp(Cmd->Keyword[i],"NChan")) 
	  sprintf(KeyValue,"%d",NChan);
	
	/* Initialize BW and CentFreq variables */
	BW=CentFreq=0.;
	for (j=0; j<NChan; j++) {
	  /* read each channel in succession */
	  fits_read_col(Fin, TDOUBLE, indx++, 1, 1, 1, NULL, 
			&ChanFreq, &anynull, &status);
	  /* Take not of first channel if StartChan keyword is requested */
	  if (!strcasecmp(Cmd->Keyword[i],"StartChan") && j==0)
	    sprintf(KeyValue,"%6.1lf",ChanFreq);
	  /* Calculate bandwidth or centre frequency using first and last
	     channels */
	  if((!strncasecmp(Cmd->Keyword[i],"BW",2) || 
	      !strncasecmp(Cmd->Keyword[i],"CFRQ",4))){
	    if (j==0)
	      BW = CentFreq = ChanFreq;
	    if (j==NChan-1) {
	      /* +4 because we are subtracting centres of channels, 
		 so this includes the edges...  */
	      BW = fabs(ChanFreq-BW)+4.;  
	      CentFreq = 0.5*(CentFreq + ChanFreq);
	    }
	  }
	  /* Print each channel separately if AllChans keyword is requestes */
	  if (!strcasecmp(Cmd->Keyword[i],"AllChan")) 
	    printf("%6.1lf    ",ChanFreq);
	}
	/* Take note of last channel if requested */
	if (!strcasecmp(Cmd->Keyword[i],"EndChan"))
	  sprintf(KeyValue,"%6.1lf",ChanFreq);
	/* Copy centre frequency and bandwidth keywords to keyword variable 
	   if requested */
	if (!strncasecmp(Cmd->Keyword[i],"CFRQ",4)) 
	  sprintf(KeyValue,"%6.1lf",CentFreq);
	if(!strncasecmp(Cmd->Keyword[i],"BW",2))
	   sprintf(KeyValue,"%6.1lf",BW);
	/* If AllChans is requested, print newline */
	if (!strcasecmp(Cmd->Keyword[i],"AllChan")) {
	  special_key=1;
	  printf("\n");
	}
	
	/* Move file pointer back to beginning of file */
	fits_movrel_hdu(Fin, -1, NULL, &status);

      /**************** End channel-based keywords: *****************/
	
      }

      /**************** DM keyword: *****************/

      else if(!strcasecmp(Cmd->Keyword[i],"DM")) {
	//	special_key=1;
	fits_movnam_hdu(Fin, ASCII_TBL, "COHDDISP", 0, &status);
	ffgnrw(Fin, &nrows, &status);
	indx = 1;
	fits_read_col(Fin, TDOUBLE, indx, 1, 1, 1, NULL, 
		      &DM, &anynull, &status);
	sprintf(KeyValue,"%8.5lf",DM);
	fits_movrel_hdu(Fin, -2, NULL, &status);

      /**************** End DM keyword: *****************/

      }
 
     /**************** NDUMPS or SCANLEN keyword: *****************/

      else if(!strcasecmp(Cmd->Keyword[i],"NDUMP") ||  // # of dumps
	      !strcasecmp(Cmd->Keyword[i],"NDUMPS") ||
	      !strcasecmp(Cmd->Keyword[i],"TDUMP") ||  // dump length in secs
	      !strcasecmp(Cmd->Keyword[i],"SCANLEN")){  // length of scans in secs

	/* Get the number of dumps.  this depensd on version of asp file code */
	fits_get_num_hdus(Fin, &NumHDU, &status);
	fits_read_key(Fin, TSTRING, "HDRVER", HdrVer, NULL, &status); 
	status=0;
	
	if(!strcasecmp(HdrVer,"Ver1.0")){
	  /* divided by 2 since half the tables are phase/period vs. freq 
	     before each data table */
	  NDump = NumHDU-3;  
	}
	else if(!strcasecmp(HdrVer,"Ver1.0.1")){
	  /* the "3" is temporary, depending on how many non-data tables we 
	     will be using */	  
	  NDump = (NumHDU-3)/2;
	}
	else{
	  printf("Do not recognize FITS file version number in header.\n");
	  printf("This header %s. Exiting...\n",HdrVer);
	  exit(3);
	}
	
	if (!strncasecmp(Cmd->Keyword[i],"NDUMP",4)) { // covers NDUMPS also...
	  sprintf(KeyValue,"%d",NDump);
	}
	// if need TDUMP or Total Scan Length (SCANLEN):
	if (!strcasecmp(Cmd->Keyword[i],"TDUMP") ||  // dump length in secs
	    !strcasecmp(Cmd->Keyword[i],"SCANLEN")) {  
	  
	  /* Get time between successive scans */

	  /* Get start time */
	  fits_read_key(Fin, TDOUBLE, "STT_SMJD", &StartTime, NULL, &status); 

	  if (NDump > 1) {

	    sprintf(aspoutstr,"ASPOUT%d",NDump-2);
	    sprintf(dumprefstr,"DUMPREF%d",NDump-2);
	    //	    printf ("\n%s   %s\n",aspoutstr,dumprefstr);
	    /* Move to first data table , get central time stamp */
	    if(!strcasecmp(HdrVer,"Ver1.0")) {
	      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT0", 0, &status);
	      fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", 
			    &(DumpMiddleSecs),NULL, &status); status = 0;
	    }
	    else if (!strcasecmp(HdrVer,"Ver1.0.1")) {
	      fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF0", 0, &status);
	      fits_read_key(Fin, TDOUBLE, "MIDSECS", &(DumpMiddleSecs), 
			    NULL, &status); status = 0;
	    }
	    TDump = DumpMiddleSecs;
	    
	    /* Move to second data table , get central time stamp */
	    sprintf(aspoutstr,"ASPOUT%d",234);
	    sprintf(dumprefstr,"DUMPREF%d",234);
	    //	    printf ("\n%s   %s\n",aspoutstr,dumprefstr);
	    if(!strcasecmp(HdrVer,"Ver1.0")) {
	      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT1", 0, &status);
	      fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", 
			    &(DumpMiddleSecs),NULL, &status); status = 0;
	    }
	    else if (!strcasecmp(HdrVer,"Ver1.0.1")) {
	      fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF1", 0, &status);
	      fits_read_key(Fin, TDOUBLE, "MIDSECS", &(DumpMiddleSecs), 
			    NULL, &status); status = 0;
	    }
	    if(DumpMiddleSecs < TDump) DumpMiddleSecs += 86400.;

	    //	    printf("First = %lf,  Second = %lf\n", TDump, DumpMiddleSecs);
	    TDump = DumpMiddleSecs - TDump;
	    
	  }
	  else {
	    /* We may have switched over to the next MJD, so take care of 
	       this if it's the case */
	    if (DumpMiddleSecs < StartTime) DumpMiddleSecs += 86400.;
	    
	    /* TDump is twice the difference between the middle time stamp 
	       and the start time */
	    TDump = 2.* ((double)floor(DumpMiddleSecs) - StartTime);
	  }
	  
	  if (!strncasecmp(Cmd->Keyword[i],"TDUMP",4)) {
	    sprintf(KeyValue,"%9.3lf",TDump);
	  }
	  if (!strcasecmp(Cmd->Keyword[i],"SCANLEN")) {
	    ScanLen = NDump*TDump;
	    sprintf(KeyValue,"%7.1lf",ScanLen);
	  }

	  /* Move back to first HDU again */

	  fits_movabs_hdu(Fin,1,&hdutype,&status);status = 0;
	}

      }


      /**************** Other keywords: *****************/
      else{
	fits_read_key(Fin, TSTRING, Cmd->Keyword[i], KeyValue, NULL, &status); 
	if(status == VALUE_UNDEFINED || !strcasecmp(KeyValue,"")) {
	  sprintf(KeyValue,"<not_found>");
	}

	/**************** End other keywords: *****************/
      }
            
      /* The only special keyword is AllChans. If it's not, do the following */
      if(!special_key){
	/* Column format */
	if(Cmd->ColumnP) {
	  printf("%s  ",KeyValue);
	  if(i==Cmd->KeywordC-1 && n_file<Cmd->InfileC-1)  printf("\n");
	}
	/* Non-column format */
	else {
	  printf("%15s:  %s\n",Cmd->Keyword[i],KeyValue);
	}
      }

    }
    status=0;

   
    fits_close_file(Fin, &status);
  }
  printf("\n");
  fflush(stdout);
  exit(0);
}
