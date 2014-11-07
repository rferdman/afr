#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPCommon.h"
#include "fitsio.h"

int ReadPSRFITSHdr(struct ASPHdr *hdr, fitsfile *Fin, struct RunVars *RunMode)
{
  int     i_chan, i_param, n_param, i_poly, i_row, retval;
  int     hdutype, colnum, anynul, status=0, tempint, dedisp;
  int     dedisp_col[16]; /* Probably not more than 16 rows to HISTORY table */
  double  ref_mjd, StartMJD, tempdouble;
  long    nrows=0;
  /* tempstr1 is for reading PSRPARAM table -- has width 128 */
  char    tempRA[32], tempDec[32], tempstr1[128], tempstr2[256];
  char    temp_par_file[128], tempo_cmd[256], tempo_cmd_0[256];
  double  hr, deg, min, sec, sgn;

  retval = -1;

  /* Move to primary header */
  if (fits_movabs_hdu(Fin, 1, &hdutype, &status)) {
    fprintf(stderr, "ERROR ReadPSRFITSHdr: Cannot move to primary HDU\n");
    return -1;
  }

/* Header version no. */
  fits_read_key(Fin, TSTRING, "HDRVER", hdr->gen.HdrVer, NULL, &status); 
  status=0; 
/* Header FITS definition type */
  fits_read_key(Fin, TSTRING, "FITSTYPE", hdr->gen.FitsType, NULL, &status); 
  status=0;

/* File creation date (YYYY-MM-DDThh:mm:ss UTC) */ 
  fits_read_key(Fin, TSTRING, "DATE", tempstr1, NULL, &status); 
  status=0;
  sscanf(tempstr1, "%sT%s", hdr->filecont.FileDate, hdr->filecont.FileUT);

/* Observer name */
  fits_read_key(Fin, TSTRING, "OBSERVER", hdr->gen.Observer, NULL, &status); 
  status=0; 
/* Project ID */
  fits_read_key(Fin, TSTRING, "PROJID", hdr->gen.ProjID, NULL, &status); 
  status=0; 
/* Observatory code */
  fits_read_key(Fin, TSTRING, "TELESCOP", hdr->obs.ObsvtyCode, NULL, &status); 
  status=0; 
  /*  if(!strncmp(telname,"GBT",3))
      sprintf(hdr->obs.ObsvtyCode, "1"); */

/* Antenna ITRF X, Y, and Z coordinates */
  fits_read_key(Fin, TDOUBLE, "ANT_X", &hdr->target.Ant_X, NULL, &status); 
  status=0; 
  fits_read_key(Fin, TDOUBLE, "ANT_Y", &hdr->target.Ant_Y, NULL, &status); 
  status=0; 
  fits_read_key(Fin, TDOUBLE, "ANT_Z", &hdr->target.Ant_Z, NULL, &status); 
  status=0; 

/* Frontend name GREG/CH etc */
  fits_read_key(Fin, TSTRING, "FRONTEND", hdr->gen.FEName, NULL, &status); 
  status=0; 

/* Beam ID for multiple-beam systems */
  fits_read_key(Fin, TSTRING, "IBEAM", hdr->gen.IBeam, NULL, &status); 
  status=0; 
  


/* Number of receiver polarisation channels */
  fits_read_key(Fin, TINT, "NRCVR", &hdr->gen.NRcvr, NULL, &status); 
  status=0;

/* Polarisation type CIRC/LIN */
  fits_read_key(Fin, TSTRING, "FD_POLN", hdr->gen.FEPol, NULL, &status); 
  status=0;
  if (!strncmp(hdr->gen.FEPol, "LIN", 3))
    strcpy(hdr->gen.FEPol, "L");
  else if (!strncmp(hdr->gen.FEPol, "CIRC", 4))
    strcpy(hdr->gen.FEPol, "C");
  else {
    fprintf(stderr,"ReadPSRFITSHdr ERROR: Unrecognized polarization basis %s\n",
	    hdr->gen.FEPol);
    return -1;
  }

/* A/B polarization swap -- +/- 1.  +1 is LIN: A=X, B=Y, CIRC:A=L, B=R */
  fits_read_key(Fin, TINT, "FD_HAND", &hdr->gen.FDHand, NULL, &status); 
  status=0; 
/* FA of E vector for equal sig in A&B */
  fits_read_key(Fin, TFLOAT, "FD_SANG", &hdr->gen.FDSang, NULL, &status); 
  status=0; 
/* Phase of A*B for injected cal */ 
  fits_read_key(Fin, TFLOAT, "FD_XYPH", &hdr->gen.FDXYph, NULL, &status); 
  status=0; 
 

/* Backend name ASP/ABPP/WAPP/PUPPI/GUPPI etc */
  fits_read_key(Fin, TSTRING, "BACKEND", hdr->gen.BEName, NULL, &status); 
  status=0; 
/* Backend Config filename */
  fits_read_key(Fin, TSTRING, "BECONFIG", hdr->gen.BEConfFile, NULL, &status); 
  status=0; 

/* 0/+1/-1 BE cross-phase: 0 unknown, +/-1 std/rev */
  fits_read_key(Fin, TINT, "BE_PHASE", &hdr->gen.BEPhase, NULL, &status); 
  status=0; 

/* 0/1 BE downconversion conjugation */
  fits_read_key(Fin, TINT, "BE_DCC", &hdr->gen.BEDownConv, NULL, &status); 
  status=0; 

/* Backend propn delay from digitiser in seconds */
  fits_read_key(Fin, TFLOAT, "BE_DELAY", &hdr->gen.BEDelay, NULL, &status); 
  status=0; 

/* On-line cycle time */
  fits_read_key(Fin, TDOUBLE, "TCYCLE", &hdr->gen.TCycle, NULL, &status); 
  status=0; 



/* Observation mode PSR/CAL/SPECTRAL etc */
  fits_read_key(Fin, TSTRING, "OBS_MODE", hdr->gen.ObsMode, NULL, &status); 
  status=0; 

/****** DATE-OBS ******/
/* Date of observation (YYYY-MM-DDThh:mm:ss UTC) */ 
  fits_read_key(Fin, TSTRING, "DATE-OBS", tempstr1, NULL, &status); 
  status=0;
  sscanf(tempstr1, "%sT%s", hdr->obs.StartDate, hdr->obs.StartUT);

/* Overall Central Sky Frequency [MHz] */
  fits_read_key(Fin, TDOUBLE, "OBSFREQ", &hdr->obs.FSkyCent, NULL, &status); 
  status=0; 
/* Bandwidth for observation [MHz] */
  fits_read_key(Fin, TDOUBLE, "OBSBW", &hdr->obs.BW, NULL, &status); 
  status=0; 
  /* Ensure BW is a positive number -- will be calculating 
     sideband later anyway*/
  hdr->obs.BW = fabs(hdr->obs.BW);
/* Number of frequency channels (original number) */
  fits_read_key(Fin, TINT, "OBSNCHAN", &hdr->obs.NChanOrig, NULL, &status); 
  status=0; 

/* DM used for online dedispaersion [cm-3 pc] */
  fits_read_key(Fin, TDOUBLE, "CHAN_DM", &hdr->obs.ChanDM, NULL, &status); 
  status=0; 

/* Name or ID for pointing centre (multibeam feeds) */
  fits_read_key(Fin, TSTRING, "PNT_ID", hdr->gen.PntID, NULL, &status); 
  status=0;   

/* Source Name */
  fits_read_key(Fin, TSTRING, "SRC_NAME", hdr->target.PSRName, NULL, &status); 
  status=0; 
/* Coordinate mode for following 4 entries  J2000/GAL/ECLIPTIC/AZEL etc. */
  fits_read_key(Fin, TSTRING, "COORD_MD", hdr->target.CoordMode, NULL, &status);
  status=0;

/* Epoch of the coordinates */
/* Read in as string and convert depending on PSRFITS version */
  fits_read_key(Fin, TFLOAT, "EQUINOX", &hdr->target.Epoch, NULL, &status); 
  status=0; 
 
/* Right Assention hh.hhhhhhhhhhhhhh */
  fits_read_key(Fin, TSTRING, "RA", tempRA, NULL, &status); 
  if (status != 0)   /* if no RA keyword, copy in STT_CRD1 keyword */
    strcpy(tempRA, hdr->target.StartCrd1); 
  sscanf(tempRA, "%lf:%lf:%lf", &hr, &min, &sec);
  hdr->target.RA = hr + min/60. + sec/3600.;
/* Declination dd.dddddddddddddd */
  fits_read_key(Fin, TSTRING, "DEC", tempDec, NULL, &status); 
  if (status != 0)   /* if no DEC keyword, copy in STT_CRD2 keyword */
    strcpy(tempDec, hdr->target.StartCrd2);
  sscanf(tempDec, "%lf:%lf:%lf", &deg, &min, &sec);
  if(deg > 0.) sgn = 1.; else sgn = -1;
  hdr->target.Dec = deg + sgn*(min/60. + sec/3600.);

/* Coordinate 1 at stop time */
  fits_read_key(Fin, TDOUBLE, "BMAJ", &hdr->target.BeamMajor, NULL, &status); 
  status=0;   
/* Beam minor axis length */
  fits_read_key(Fin, TDOUBLE, "BMIN", &hdr->target.BeamMinor, NULL, &status); 
  status=0;   
/* Beam position angle */
  fits_read_key(Fin, TDOUBLE, "BPA", &hdr->target.BeamPA, NULL, &status); 
  status=0;   

/* coordinate 1 at start time */
  fits_read_key(Fin, TSTRING, "STT_CRD1",hdr->target.StartCrd1, NULL, &status);
  status=0; 
/* coordinate 2 at start time */
  fits_read_key(Fin, TSTRING, "STT_CRD2",hdr->target.StartCrd2, NULL, &status);
  status=0; 
/* coordinate 1 at stop time */
  fits_read_key(Fin, TSTRING, "STP_CRD1", hdr->target.StopCrd1, NULL, &status);
  status=0; 
/* coordinate 2 at stop time */
  fits_read_key(Fin, TSTRING, "STP_CRD2", hdr->target.StopCrd2, NULL, &status);
  status=0; 

  
/* Track mode TRACK/SLEW/DRIFT/SCANGO/SCANLAT etc. */
  fits_read_key(Fin, TSTRING, "TRK_MODE", hdr->target.TrackMode, NULL, &status);
  status=0; 
/* Requested observation length (s) */
  fits_read_key(Fin, TDOUBLE, "SCANLEN", &hdr->obs.ObsLength, NULL, &status); 
  status=0; 

/* Feed track mode -- FA/CPA/SPA/TPA */
  fits_read_key(Fin, TSTRING, "FD_MODE", hdr->gen.FeedMode, NULL, &status);
  status=0; 
/* Feed/Position angle requested */
  fits_read_key(Fin, TDOUBLE, "FA_REQ", &hdr->gen.FAReq, NULL, &status);
  status=0; 

/* Cal mode (OFF, SYNC, EXT1, EXT2) */
  fits_read_key(Fin, TSTRING, "CAL_MODE", hdr->gen.CalMode, NULL, &status);
  status=0; 
/* Cal modulation frequency in Hz */  
  fits_read_key(Fin, TDOUBLE, "CAL_FREQ", &hdr->obs.CalFreq, NULL, &status);
  status=0; 
/* Cal duty cycle */  
  fits_read_key(Fin, TDOUBLE, "CAL_DCYCLE", &hdr->obs.CalDutyCycle, NULL, &status);
  status=0; 
/* Cal phase (wrt start time) */
  fits_read_key(Fin, TDOUBLE, "CAL_PHS", &hdr->obs.CalPhase, NULL, &status);
  status=0; 
/* Number of states in cal pulse */
  fits_read_key(Fin, TDOUBLE, "CAL_NPHS", &hdr->obs.CalNStates, NULL, &status);
  status=0; 
    

/* Integer MJD of starting sample */
  fits_read_key(Fin, TINT,    "STT_IMJD", &hdr->obs.IMJDStart, NULL, &status); 
  status=0;
/* Start Time (INT) second  past 0.0h */
  fits_read_key(Fin, TINT, "STT_SMJD", &hdr->obs.StartTime, NULL, &status); 
  status=0; 
 
/* Clock offset (sec) */
  fits_read_key(Fin, TDOUBLE, "STT_OFFS", &hdr->obs.ClockOffset, NULL, &status);
  status=0; 
/* Start LST */
  fits_read_key(Fin, TDOUBLE, "STT_LST", &hdr->target.StartLST, NULL, &status); 
  status=0; 

/****  Get DM/RM infor ****/
  /* Don't bother to get this DM/RM info if this is a cal file */
  if (strncmp(hdr->gen.ObsMode, "CAL", 3) &&
      strncmp(hdr->gen.ObsMode, "FOF", 3) &&
      /* Not sure what FON is, but got this bringing up errors in
	 this if statement, so ignore if ObsMode is "FON"... */
      strncmp(hdr->gen.ObsMode, "FON", 3) &&
      strncmp(hdr->gen.ObsMode, "SEARCH", 6) ) {
      
      /* If this is GUPPI data, we need to read in DM/RM from the ephemeris table */
    //    if (!strncmp(hdr->gen.BEName, "GUPPI", 5)){
      /* Go to PSREPHEM table -- if nonexistent, see if there is a PSRPARAM
	 table to extract DM */
      
      if(fits_movnam_hdu(Fin, BINARY_TBL, "PSREPHEM", 0, &status)) {
	/* Now try PSRPARAM table */
	status=0;
	if(fits_movnam_hdu(Fin, BINARY_TBL, "PSRPARAM", 0, &status)) {
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: neither PSREPHEM nor PSRPARAM tables exist!\n");
	  return -1;
	}
	else {
	  /* Need to read row by row until DM is seen -- and 
	     make sure it's not DM1, etc.*/
	  fits_get_num_rows (Fin, &nrows, &status);
	  if(!nrows){
	    fprintf(stderr, "ERROR ReadPSRFITSHdr: No rows in PSRPARAM table.\n");
	    return -1;
	  }

	  /* First get number of rows */
	  fits_read_key(Fin, TINT, "NAXIS2", &n_param, NULL, &status); 
	  if(nrows!=n_param){
	    fprintf(stderr, "ERROR ReadPSRFITSHdr: Problem with row numbers in PSRPARAM table.\n");
	    return -1;
	  }
	  /*  will need to read in DM and RM from here as a table column... */
	  if(fits_get_colnum(Fin, CASEINSEN, "PARAM", &colnum, &status)){
	    fprintf(stderr, "ERROR ReadPSRFITSHdr: Problem reading PSRPARAM table.\n");
	    return -1;
	  }

	  char* temp = &(tempstr1[0]);
	  for (i_param=1; i_param <= n_param; i_param++) {
	    //#if 0
	    if(fits_read_col(Fin, TSTRING, colnum, i_param, 1, 1, 
			     0, &temp, &anynul, &status)) {
	      fprintf(stderr, "ERROR ReadPSRFITSHdr: Problem reading PSRPARAM table.\n");
	      return -1;
	    }
	    //#endif


	    /* Parse string */
	    printf("%s\n",tempstr1);
	    if(!strncmp(&tempstr1[0], "DM ", 3)){
	      sscanf(tempstr1, "%s %lfD%3d", tempstr2, &tempdouble, &tempint);
	      sprintf(tempstr2,"%lfE%d   ", tempdouble, tempint);
	      sscanf(tempstr2, "%lf", &hdr->obs.DM);
	      /* Need to set RM to zero in this case */
	      hdr->obs.RM = 0.;
	      break;
	    }
	    if(i_param==n_param-1){
	      fprintf(stderr, "ERROR ReadPSRFITSHdr:  Could not find DM in PSRPARAM table.\n");
	      return -1;
	    }	    
	  }

	}
      }
      else{
	
	/*  will need to read in DM and RM from here as a table column... */
	if(fits_get_colnum(Fin, CASEINSEN, "DM", &colnum, &status)){
	  fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read DM column.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, 1, NULL, &hdr->obs.DM, &anynul, &status)) {
	  fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read DM value.\n");
	}
	
	if(fits_get_colnum(Fin, CASEINSEN, "RM", &colnum, &status)){
	  fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read RM column.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, 1, NULL, &hdr->obs.RM, &anynul, &status)) {
	  fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read RM value.\n");
	}
      }
      // }
  }
  
  /* Will need to find out whether the PSRFITS profile data are already 
     dedispersed */

  /* Initialize header variable to zero, i.e. data has not already been 
     dedispersed */
  hdr->redn.Dedisp = 0;

  /* Go to HISTORY table */
  if(fits_movnam_hdu(Fin, BINARY_TBL, "HISTORY", 0, &status)) {
    fprintf(stderr, "ReadASPHdr WARNING: HISTORY table does not exist!\n");
    // return -1;
	status = 0;
  } 
  else{
	  /* Get number of rows */
	  fits_get_num_rows(Fin, &nrows, &status);
	  if(nrows > 0){
		  /* Read DEDISP table, and get the last row's value */
		  if(fits_get_colnum(Fin, CASEINSEN, "DEDISP", &colnum, &status)){
			  fprintf(stderr, "ReadPSRFITSHdr ERROR: Problem reading DEDISP column ");
			  fprintf(stderr, "in HISTORY table.\n");
			  return -1;
		  }
		  if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, 1, NULL, dedisp_col, &anynul, &status)) {
			  fprintf(stderr, "ReadPSRFITSHdr ERROR: Could not read DEDISP column ");
			  fprintf(stderr, "in HISTORY table.\n");
			  /* In this case, assume no dedispersion was done by default */
			  hdr->redn.Dedisp = 0; 
		  }
		  /* Sum DEDISP column values */
		  dedisp=0;
		  for (i_row=0; i_row<nrows; i_row++){
			  dedisp += dedisp_col[i_row];
		  }
		  if(dedisp > 0){ /* i.e. data HAS been dedispersed already */
			  hdr->redn.Dedisp = 1;
			  /* Turn off dedispersion, overriding possible command line input */
			  printf("Data is already dedispersed.\n");
			  if(RunMode->Dedisp) printf("Turning off dedispersion request.\n");
			  RunMode->Dedisp = 0;
		  }
	  }
	  /* If there are no rows or is no HISTORY table, assume no dedispersion 
	  was done -- header variable remains at zero */
  }


  /* Dispersion Measure */
    //  fits_read_key(Fin, TDOUBLE, "DM", &hdr->obs.DM, NULL, &status); status=0; 
  /* Rotation Measure */
  //fits_read_key(Fin, TDOUBLE, "RM", &hdr->obs.RM, NULL, &status); status=0; 

  /* Will need to get data info from different headers and tables,
     depending on version of PSRFITS... */

  if (!strncmp(hdr->gen.HdrVer, "1.26", 4)) {

    /* Go to HISTORY table */
    if(fits_movnam_hdu(Fin, BINARY_TBL, "HISTORY", 0, &status)) {
      fprintf(stderr, "ERROR ReadASPHdr: HISTORY table does not exist!\n");
      return -1;
    } 
    
   /* Number of profile bins */
    if(fits_get_colnum(Fin, CASEINSEN, "NBIN", &colnum, &status)){
      fprintf(stderr,"ERROR ReadPSRFITSHdr: Could not read NBIN column.\n");
      return -1;
    }
    if(fits_read_col(Fin, TINT, colnum, 1, 1, 1, NULL, 
		     &hdr->redn.RNBinTimeDump, &anynul, &status)) {
      fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read NBIN value.\n");
    }
   
    /* Number of frequency channels (current file) */
    if(fits_get_colnum(Fin, CASEINSEN, "NCHAN", &colnum, &status)){
      fprintf(stderr,"ERROR ReadPSRFITSHdr: Could not read NCHAN column.\n");
      return -1;
    }
    if(fits_read_col(Fin, TINT, colnum, 1, 1, 1, NULL, 
		     &hdr->obs.NChan, &anynul, &status)) {
      fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read NCHAN value.\n");
      return -1;
    }

   /* Channel/sub-band width (current file) */
    if(fits_get_colnum(Fin, CASEINSEN, "CHAN_BW", &colnum, &status)){
      fprintf(stderr,"ERROR ReadPSRFITSHdr: Could not read CHAN_BW column.\n");
      return -1;
    }
    if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, 1, NULL, 
		     &hdr->obs.ChanWidth, &anynul, &status)) {
      fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read CHAN_BW value.\n");
      return -1;
    }
    /* Ensure ChanWidth is positive -- will calculate sideband later anyway */
    hdr->obs.ChanWidth = fabs(hdr->obs.ChanWidth);

    /* Number of polarizations */
    if(fits_get_colnum(Fin, CASEINSEN, "NPOL", &colnum, &status)){
      fprintf(stderr,"ERROR ReadPSRFITSHdr: Could not read NPOL column.\n");
      return -1;
    }
    if(fits_read_col(Fin, TINT, colnum, 1, 1, 1, NULL, 
		     &hdr->obs.NPoln, &anynul, &status)) {
      fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read NPOL value.\n");
      return -1;
    }

    /* Polarization identifier */
    if(fits_get_colnum(Fin, CASEINSEN, "POL_TYPE", &colnum, &status)){
      fprintf(stderr,"ERROR ReadPSRFITSHdr: Could not read POL_TYPE column.\n");
      return -1;
    }
    //int i_char;
    char* temppoltype = &(hdr->obs.PolnType[0]);
    fits_read_col(Fin, TSTRING, colnum, 1, 1, 1, NULL, 
		  &temppoltype, &anynul, &status);
    if(status!=0){
      fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read POL_TYPE value.\n");
      fits_report_error(stderr, status);
      return -1;
    }
    
    /* Go to SUBINT table */
    if(fits_movnam_hdu(Fin, BINARY_TBL, "SUBINT", 0, &status)) {
      fprintf(stderr, "ReadASPHdr ERROR: SUBINT table does not exist!\n");
      return -1;
    } 
  }
  else{

  /* Go to SUBINT table */
    if(fits_movnam_hdu(Fin, BINARY_TBL, "SUBINT", 0, &status)) {
      fprintf(stderr, "ReadASPHdr ERROR: SUBINT table does not exist!\n");
      return -1;
    } 
    /* First, if not GUPPI data, read in DM and RM from this table's header */
    if(strncmp(hdr->gen.BEName, "GUPPI", 5)) {
      fits_read_key(Fin, TDOUBLE, "DM", &hdr->obs.DM, NULL, &status); 
      status=0; 
      /* Rotation Measure */
      fits_read_key(Fin, TDOUBLE, "RM", &hdr->obs.RM, NULL, &status); 
      status=0; 
    }
    
    /***** TBIN -- time per bin or sample *****/
    
    /* Number of profile bins */
    fits_read_key(Fin, TINT, "NBIN", &hdr->redn.RNBinTimeDump, NULL, &status);
    status=0;
    /* Number of frequency channels (current file) */
    fits_read_key(Fin, TINT, "NCHAN", &hdr->obs.NChan, NULL, &status); 
    status=0; 
    /* Channel/sub-band width (current file) */
    fits_read_key(Fin, TDOUBLE, "CHAN_BW", &hdr->obs.ChanWidth, NULL, &status);
    status=0;
    /* Ensure ChanWidth is positive -- will calculate sideband later anyway */
    hdr->obs.ChanWidth = fabs(hdr->obs.ChanWidth);
    /* Number of polarizations */
    fits_read_key(Fin, TINT, "NPOL", &hdr->obs.NPoln, NULL, &status); 
    status=0; 
    /* Polarization identifier */
    fits_read_key(Fin, TSTRING, "POL_TYPE", hdr->obs.PolnType, NULL, &status); 
    status=0; 
    /* Read in channel values */
    /* NOTE that there is a separate DAT_FREQ table for each row...  
       hopefully will be the same for all integrations!  
       Will put in a safeguard at some point */
    if(fits_read_key(Fin, TINT, "NSUBOFFS", &hdr->obs.NSubOffs, NULL, &status)){
      if(!strncmp(hdr->gen.ObsMode, "CAL", 3)){
		  fprintf(stderr, "ReadPSRFITSHdr WARNING:  Could not get NSUBOFFS from SUBINT extension.\n");
		  printf("Proceeding with CAL file reading.\n");	
		  status=0;  /* Otherwise next fits routine will crash */
	
      } 
      else {
		  fprintf(stderr, "ReadPSRFITSHdr WARNING:  Could not read NSUBOFFS keyword from SUBINT extension.\n");
		  printf("Proceeding with reading of data.  Time stamps likely not correct.\n");	
		  hdr->obs.NSubOffs=0;
		  //sleep(2);
		  status=0;  /* Otherwise next fits routine will crash */
		  //  fits_report_error(stderr, status); /* print any error message */
	
      }
    }
    
  }
  
  /* Number of integrations */
  fits_read_key(Fin, TINT, "NAXIS2", &hdr->redn.RNTimeDumps, NULL, &status); 
  status=0;
  
/* Length of integration */
  if(fits_get_colnum (Fin, CASEINSEN, "TSUBINT", &colnum, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: No TSUBINT in FITS file?\n");
    return -1;
  }
  if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, 1, NULL, &hdr->redn.TDump, &anynul, &status)){
    fprintf(stderr, "ReadPSRFITSData ERROR: Unable to read TSUBINT...\n");
    return -1;
  }

  if(fits_get_colnum(Fin, CASEINSEN, "DAT_FREQ", &colnum, &status)){
    fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read DAT_FREQ column.\n");
    return -1;
  }
  if(fits_read_col(Fin, TDOUBLE, colnum, 1, 1, hdr->obs.NChan, NULL, &hdr->obs.ChanFreq[0], &anynul, &status)) {
    fprintf(stderr, "ERROR ReadPSRFITSHdr: Could not read DAT_FREQ value.\n");
    return -1;
  }
  

  /* Set BW -- should work for all backend versions of PSRFITS */
  hdr->obs.BW = hdr->obs.ChanWidth * ((double)hdr->obs.NChan);
  /* Set average Frequency here -- should work for all versions of PSRFITS */
  //  hdr->obs.FSkyCent = 0.5*(hdr->obs.ChanFreq[0] + hdr->obs.ChanFreq[hdr->obs.NChan-1]);
  hdr->obs.ChanFreqMean = 0.5*(hdr->obs.ChanFreq[0] + hdr->obs.ChanFreq[hdr->obs.NChan-1]);




  /* Now read in polycos if needed, and if POLYCO table exists, to get 
     individual profile phase information.  This can also be used to shift 
     profiles later on using new ephemeris or set of polycos, if given 
     on command line */

  /* First, allocate Header Polyco array: */
  if(RunMode->ForcePoly)
    //  if(!strcmp(hdr->gen.BEName, "DFB") || !strcmp(hdr->gen.BEName, "Jod"))
    /* Allocate Polyco structure to number of channels */
    /* For now, only needed if user has specified -forcepoly on command line */
    hdr->redn.Polycos = (struct Polyco *)malloc(hdr->obs.NChan*MAX_PC_SETS*
						sizeof(struct Polyco));  
  else
    /* Allocate Polyco structure to only deal with one channel */
    hdr->redn.Polycos = (struct Polyco *)malloc(MAX_PC_SETS*
						sizeof(struct Polyco));    
  
  /* First, (try to) move to POLYCO extension */

  /* Only do for observation mode PSR */
  if(!strcmp(hdr->gen.ObsMode, "PSR")) {
    if(fits_movnam_hdu(Fin, BINARY_TBL, "POLYCO", 0, &status)) {
      /* POLYCO table not found.  Set RunMode variable to 0 */
      fprintf(stderr, "ReadPSRFITSHdr WARNING: POLYCO table does not exist\n");
      RunMode->PolyTable = 0;
      /* Reset status = 0 so that further fitsio routines will work */
      status = 0;
      /* Assign reference frequency to be centre frequency for the purposes of 
	 initial alignment */
      hdr->redn.NPoly = 0;
      hdr->redn.Polycos[0].FSkyRef = hdr->obs.FSkyCent;
    }
    else {
      /* POLYCO table found -- Set RunMode variable to 1 */
      RunMode->PolyTable = 1;
    }
    
    /* IF POLYCO table was found, read it in and assign values to Polyco 
       structure, including coefficients */
    if(RunMode->PolyTable){
      /* Read in number of rows, ie. number of polyco sets */
      if(fits_read_key(Fin, TINT, "NAXIS2", &hdr->redn.NPoly, NULL, &status)) {
	fprintf(stderr, "ReadPSRFITSHdr ERROR:  Could not get NAXIS2 from POLYCO extension.");
	return -1;
      }
      /* Now read in actual information for each polyco set in data file */
      /* NSPAN -- Valid span of polyco set in minutes */
      for(i_poly=0; i_poly<hdr->redn.NPoly; i_poly++){
	if(fits_get_colnum (Fin, CASEINSEN, "NSPAN", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No NSPAN in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TINT, colnum, 1+i_poly, 1, 1, NULL, &hdr->redn.Polycos[i_poly].NMinutes, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read NSPAN...\n");
	  return -1;
	}
	
	/* NCOEFF -- Number of coefficients*/
	if(fits_get_colnum (Fin, CASEINSEN, "NCOEF", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No NCOEF in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TINT, colnum, 1+i_poly, 1, 1, NULL, 
			 &hdr->redn.Polycos[i_poly].NCoeff, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read NCOEF...\n");
	  return -1;
	}
	
	/* REF_MJD -- Reference MJD for this polyco set */
	if(fits_get_colnum (Fin, CASEINSEN, "REF_MJD", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No REF_MJD in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, &ref_mjd, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read REF_MJD...\n");
	  return -1;
	}
	hdr->redn.Polycos[i_poly].MjdMidInt = floor(ref_mjd);
	hdr->redn.Polycos[i_poly].MjdMidFrac = ref_mjd - floor(ref_mjd);
	
	if(RunMode->Verbose){
	  printf("ref_mjd:                %lf\n", ref_mjd);
	  printf("Polycos:  MjdMidInt:    %lf\n",  
		 hdr->redn.Polycos[i_poly].MjdMidInt);
	  printf("Polycos:  MjdMidFrac:   %lf\n\n",  
		 hdr->redn.Polycos[i_poly].MjdMidFrac);
	}
	
	/* REF_FREQ -- Reference observing frequency for this polyco set */
	if(fits_get_colnum (Fin, CASEINSEN, "REF_FREQ", &colnum, &status)){
	  fprintf(stderr,"ReadPSRFITSHdr ERROR: No REF_FREQ in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, 
			 &hdr->redn.Polycos[i_poly].FSkyRef, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read REF_FREQ...\n");
	  return -1;
	}
	
	/* REF_PHS -- Reference phase for this polyco set */
	if(fits_get_colnum (Fin, CASEINSEN, "REF_PHS", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No REF_PHS in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, 
			 &hdr->redn.Polycos[i_poly].PhRotRef, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read REF_PHS...\n");
	  return -1;
	}
	
	/* REF_F0 -- Reference rotation frequency for this polyco set */
	if(fits_get_colnum (Fin, CASEINSEN, "REF_F0", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No REF_F0 in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 1, NULL, 
			 &hdr->redn.Polycos[i_poly].FRotRef, &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read REF_F0...\n");
	  return -1;
	}
	
	/* COEFF -- Actual polynomial coefficients for this set */
	if(fits_get_colnum (Fin, CASEINSEN, "COEFF", &colnum, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: No COEFF in POLYCO table.\n");
	  return -1;
	}
	if(fits_read_col(Fin, TDOUBLE, colnum, 1+i_poly, 1, 
			 hdr->redn.Polycos[i_poly].NCoeff, NULL, 
			 &hdr->redn.Polycos[i_poly].Coeff[0], &anynul, &status)){
	  fprintf(stderr, "ReadPSRFITSHdr ERROR: Unable to read REF_F0...\n");
	  return -1;
	}
      }
  

      
    }
    else{
      /* If POLYCO table was not found, read in pulsar ephemeris if found. */
      printf("POLYCO table was not found.\n");
      printf("Will now create polycos from ephemeris contained in FITS file...\n");
      /* Try to find par file extension -- either PSRPARAM or PSREPHEM -- and 
	 write it to a temporary par file from which to calculate polycos 
	 corresponding to this observation */
      sprintf(&temp_par_file[0], "temp.par");
      if (WrtPSRFITSPar(Fin, temp_par_file) < 0){
	/* If pulsar ephemeris not found, then do not do pulse realignment 
	   at all. */
	fprintf(stderr, "Problem getting par file from input FITS File %s.",
		RunMode->Infile);
	fprintf(stderr, "Continuing, with no pulse phase information.\n");
	//	return -1;
	RunMode->EphTable = 0;
      }
      else{
	RunMode->EphTable = 1;
	/* Run tempo on output par file to get polycos.  Default is to calculate
	   a single polyco set for the entire band for GUPPI/PUPPI, and ROACH 
	   data. For DFB data, will calculate multi-channel polycos, since they 
	   seem to have used tempo2 predictors. Assume tempo1 par files have 
	   been used until otherwise needed. */
	StartMJD = (double)hdr->obs.IMJDStart + 
	  ((double)hdr->obs.StartTime)/86400.;
	sprintf(tempo_cmd_0,
		"tempo -f %s -Zpsr=%s -Ztobsh=%lf -Zstart=%lf -Zspan=%d -Zsite=%s",
	     /* Channel frequency must be 3 decimal places for tempo1 predictors */
		temp_par_file, hdr->target.PSRName, 
	     /* allow an extra 2 hours of polycos */
		2. + (hdr->redn.TDump*((double)hdr->redn.RNTimeDumps)/3600.), 
	     /* Start polyco set to be 1 hour early, (and so end 1 hour later) */
		hdr->obs.IMJDStart + ((double)hdr->obs.StartTime + 
			          floor((double)hdr->obs.NSubOffs*hdr->redn.TDump)
				      - 3600.)/86400.,
		/* default at 15 minutes valid span */
		1800, hdr->obs.ObsvtyCode);
	/* Now add on frequency/ies, depending on whether user wants to force 
	 multi-channel polycos */
	if(RunMode->ForcePoly){
	  //	if(!strcmp(hdr->gen.BEName, "DFB") || !strcmp(hdr->gen.BEName, "Jod")){
	  for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
	    sprintf(tempo_cmd, "%s -Zfreq=%lf", 
		    tempo_cmd_0, hdr->obs.ChanFreq[i_chan]);
	    system(tempo_cmd);
	    /* Dump current channel's polycos to a larger file for safe keeping */
	    system("cat polyco.dat >> temp_polyco.dat");
	    /* Read polyco file and get parameter values */
	    //	for(i_chan_in=0;i_chan_in<InHdr.obs.NChan;i_chan_in++){
	    if((hdr->redn.NPoly=GetPoly("polyco.dat", hdr->target.PSRName, 
					&hdr->redn.Polycos[i_chan*MAX_PC_SETS], 
					hdr->obs.ChanFreq[i_chan], 
					StartMJD)) < 1) {
	      printf("ReadPSRFITSHdr WARNING: Could not find polycos for ");
	      printf("original data profiles at %lf MHz observing frequency.\n",
		     hdr->obs.ChanFreq[i_chan]);
	    }
	  }
	  remove("polyco.dat");
	}
	else{
	  sprintf(tempo_cmd, "%s -Zfreq=%lf", 
		  tempo_cmd_0, hdr->obs.FSkyCent);
	  system(tempo_cmd);
	  if((hdr->redn.NPoly=GetPoly("polyco.dat", hdr->target.PSRName, 
				      &hdr->redn.Polycos[0], 
				      hdr->obs.FSkyCent, 
				      StartMJD)) < 1) {
	    printf("ReadPSRFITSHdr WARNING: Could not find polycos for ");
	    printf("original data profiles at %lf MHz observing frequency.\n",
		   hdr->obs.ChanFreq[i_chan]);
	  }
	  /* Copy over polyco.dat file to a new file for safe keeping */
	  system("cp polyco.dat temp_polyco.dat");
	  remove("polyco.dat");
	}            	
	
	//	}
	/* If we made it out of the above for-loop alive, then we have succeeded
	   in generating and reading in all the polycos we want. */
	printf("Polycos succesfully found and read for original profile data.\n");
	
	/* Clean up -- remove temporary par and polyco files */
	printf("Removing temporary ephemeris and polyco files...\n");
	//remove("temp.par");
	// remove("temp_polyco.dat");	

      }
    }
  }





#if 0

  /* Little hard-coded exception for Nancay data taken before MJD 53686.05 
     -- Those channels are labelled 1 MHz too high */
  if (!strcmp(hdr->obs.ObsvtyCode, "f") && 
     ((double)hdr->obs.IMJDStart + ((double)hdr->obs.StartTime/86400.0)) 
     <= 53686.05) {
    printf("Adjusting Nancay data frequency labels by -1.0 MHz...\n");
    for (i=0; i<hdr->obs.NChan; i++) {
      if (hdr->obs.ChanFreq[i] >= 1350.0 && hdr->obs.ChanFreq[i] <= 1450.0) 
	hdr->obs.ChanFreq[i] -= 1.0;
    }
  }

#endif

  retval = 0;

  return retval;
}
