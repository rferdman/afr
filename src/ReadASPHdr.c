#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPHeader.h"
#include "fitsio.h"

int ReadASPHdr(struct      ASPHdr *hdr,
		fitsfile   *Fin)
{
  int     i, hdutype, status, knul=0, anynull, indx, retval;
  double  dnul=0.0;
  float   fnul=0;
  long    nrows=0;
  long    NPtsProf=0;
  char    tblname[40]; //, *ttype[10], *tform[10], *tunit[10];

  int     NumHDU;
  double  DumpMiddleSecs;
  
  
  retval = -1;

  fits_read_key(Fin, TSTRING, "SCANNAME", hdr->gen.ScanName, NULL, &status); 
  status=0;
/* Ver of software used for analysis */
  fits_read_key(Fin, TSTRING, "SOFT_VER", hdr->gen.SoftVer,  NULL, &status); 
  status=0; 
/* Computer platform used */  
  fits_read_key(Fin, TSTRING, "PLATFORM", hdr->gen.Platform, NULL, &status); 
  status=0;  
/* Comments by Operator/Observer */
  fits_read_key(Fin, TSTRING, "COMMOPER", hdr->gen.CommentOper, NULL, &status); 
  status=0;  
/* Header version no. */
  fits_read_key(Fin, TSTRING, "HDRVER", hdr->gen.HdrVer, NULL, &status); 
  status=0; 
/* Observer name */
  fits_read_key(Fin, TSTRING, "OBSERVER", hdr->gen.Observer, NULL, &status); 
  status=0; 
/* Project ID */
  fits_read_key(Fin, TSTRING, "PROJID", hdr->gen.ProjID, NULL, &status); 
  status=0; 
/* Observatory code */
  fits_read_key(Fin, TSTRING, "OBSVTY", hdr->obs.ObsvtyCode, NULL, &status); 
  status=0; 
/* Frontend name GREG/CH etc */
  fits_read_key(Fin, TSTRING, "FRONTEND", hdr->gen.FEName, NULL, &status); 
  status=0; 
/* Polarisation type CIRC/LIN */
  fits_read_key(Fin, TSTRING, "FD_POLN", hdr->gen.FEPol, NULL, &status); 
  status=0; 
/* Backend name ASP/ABPP/WAPP etc */
  fits_read_key(Fin, TSTRING, "BACKEND", hdr->gen.BEName, NULL, &status); 
  status=0; 
/* Backend Config filename */
  fits_read_key(Fin, TSTRING, "BECONFIG", hdr->gen.BEConfFile, NULL, &status); 
  status=0; 
/* Observation mode PSR/CAL/SPECTRAL etc */
  fits_read_key(Fin, TSTRING, "OBS_MODE", hdr->gen.ObsMode, NULL, &status); 
  status=0; 
/* Source Name */
  fits_read_key(Fin, TSTRING, "SRC_NAME", hdr->target.PSRName, NULL, &status); 
  status=0; 
/* Right Assention hh.hhhhhhhhhhhhhh */
  fits_read_key(Fin, TDOUBLE, "RA", &hdr->target.RA, NULL, &status); 
  status=0; 
/* Declination dd.dddddddddddddd */
  fits_read_key(Fin, TDOUBLE, "DEC", &hdr->target.Dec, NULL, &status); 
  status=0; 
/* Epoch of the coordinates */
  fits_read_key(Fin, TFLOAT, "EPOCH", &hdr->target.Epoch, NULL, &status); 
  status=0; 
/* Coordinate mode for following 4 entries  J2000/GAL/ECLIPTIC/AZEL etc. */
  fits_read_key(Fin, TSTRING, "COORD_MD", hdr->target.CoordMode, NULL, &status);
  status=0; 
/* coordinate 1 at start time */
  fits_read_key(Fin, TDOUBLE, "STT_CRD1",&hdr->target.StartCrd1, NULL, &status);
  status=0; 
/* coordinate 2 at start time */
  fits_read_key(Fin, TDOUBLE, "STT_CRD2",&hdr->target.StartCrd2, NULL, &status);
  status=0; 
/* coordinate 1 at stop time */
  fits_read_key(Fin, TDOUBLE, "STP_CRD1", &hdr->target.StopCrd1, NULL, &status);
  status=0; 
/* coordinate 2 at stop time */
  fits_read_key(Fin, TDOUBLE, "STP_CRD2", &hdr->target.StopCrd2, NULL, &status);
  status=0; 
/* Track mode TRACK/SLEW/DRIFT etc. */
  fits_read_key(Fin, TSTRING, "TRK_MODE", hdr->target.TrackMode, NULL, &status);
  status=0; 
/* Observation length (s) */
  fits_read_key(Fin, TDOUBLE, "SCANLEN", &hdr->obs.ObsLength, NULL, &status); 
  status=0; 
/* Start Date/Month/Year */
  fits_read_key(Fin, TSTRING, "STT_DATE", hdr->obs.StartDate, NULL, &status); 
  status=0; 
/* Start Time in UT */
  fits_read_key(Fin, TSTRING, "STT_TIME", hdr->obs.StartUT, NULL, &status); 
  status=0; 
/* Integer MJD of starting sample */
  fits_read_key(Fin, TINT,    "STT_IMJD", &hdr->obs.IMJDStart, NULL, &status); 
  status=0; 
/* Start Time (INT) past 0.0h */
  fits_read_key(Fin, TINT, "STT_SMJD", &hdr->obs.StartTime, NULL, &status); 
  status=0; 
/* Fractional second of start time */
  fits_read_key(Fin, TDOUBLE, "STT_FRAC", &hdr->obs.StartFSec, NULL, &status); 
  status=0; 
/* Clock offset (sec) */
  fits_read_key(Fin, TDOUBLE, "STT_OFFS", &hdr->obs.ClockOffset, NULL, &status);
  status=0; 
/* Start LST */
  fits_read_key(Fin, TDOUBLE, "STT_LST", &hdr->target.StartLST, NULL, &status); 
  status=0; 
  fits_read_key(Fin, TFLOAT, "ION_RM", &hdr->obs.IonRM, NULL, &status); 
  status=0; 
/* Sampling interval (sec) */
  fits_read_key(Fin, TDOUBLE, "SAMP_INT",&hdr->obs.SampInterval, NULL, &status);
  status=0; 
/* Output data type String (IQUVLT etc) */
  fits_read_key(Fin, TSTRING, "OP_STRNG", hdr->obs.OPString, NULL, &status); 
  status=0; 
/* Output scale 1-Arbit; 2-Jy; 3-K;  */
  fits_read_key(Fin, TINT, "OP_SCALE", &hdr->obs.OPScale, NULL, &status); 
  status=0; 
/* Dynamic DC subtracted? 0/1 */
  fits_read_key(Fin, TINT, "DYN_DC", &hdr->obs.DynDC, NULL, &status); 
  status=0; 
/* Dynamic DC correction time scale (s) */
  fits_read_key(Fin, TFLOAT, "DYN_TIME", &hdr->obs.DynTime, NULL, &status); 
  status=0; 
/* No. of bits per sample */
  fits_read_key(Fin, TINT, "NBITS", &hdr->obs.NBitPerSamp, NULL, &status); 
  status=0; 
/* Interference corrected? 0/1 */
  fits_read_key(Fin, TINT, "INTERF", &hdr->obs.InterfString, NULL, &status); 
  status=0; 
/* Zap filename */
  fits_read_key(Fin, TSTRING, "ZAPFILE", hdr->obs.ZapFileName, NULL, &status); 
  status=0; 
/* Number of bins per profile */
  fits_read_key(Fin, TINT, "NPTSPROF", &hdr->redn.RNBinTimeDump, NULL, &status);
  status=0;  
/* Number of bins per profile */
  fits_read_key(Fin, TINT, "NDUMPS", &hdr->redn.RNTimeDumps, NULL, &status); 
  status=0; 
/* Overall Central Sky Frequency */
  fits_read_key(Fin, TDOUBLE, "FSKYCENT", &hdr->obs.FSkyCent, NULL, &status); 
  status=0; 

  /*  for (i = 0; i < 10; i++)
    {
      ttype[i] = (char *) malloc(20);
      tform[i] = (char *) malloc(20);
      tunit[i] = (char *) malloc(20);
      }  */


  strcpy(tblname, "BECONFIG");

  fits_movnam_hdu(Fin, ASCII_TBL, "BECONFIG", 0, &status);
  fits_get_num_rows(Fin, &nrows, &status);

  indx = 1;
  fits_read_col(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.NChan, &anynull, &status);
  indx++;
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.ChanFreq[i], 
	  &anynull, &status);
    indx++;  
/*     printf("CHANFREQ %d = %f\n",i,hdr->obs.ChanFreq[i]);fflush(stdout); */
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.ChanWidth, &anynull, 
	  &status);
    indx++;
  }
  /*
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TFLOAT, indx, 1, 1, 1, &fnul, &hdr->obs.ChanTSys[i], &anynull, 
	  &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TFLOAT, indx, 1, 1, 1, &fnul, &hdr->obs.ChanTCal[i], &anynull, 
	  &status);
    indx++;
  }
  fits_read_col(Fin, TINT, indx, 1, 1, 1, &fnul, &hdr->asp.SplitFac, &anynull, 
  &status); */

  fits_movnam_hdu(Fin, ASCII_TBL, "COHDDISP", 0, &status);

  indx = 1;
  fits_read_col(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.DM, &anynull, &status); 
  indx++;
/*   printf("DM BEFORE = %f\n",hdr->obs.DM); */
  fits_read_col(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.DMMethod, &anynull, &status); 
  indx++;
  fits_read_col(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.RM, &anynull, &status); 
  indx++;
  fits_read_col(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.RMMethod, &anynull, &status); 
  indx++;

  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.ChanChirpLen[i], &anynull,
	  &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_read_col(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.ChanOverlap[i], &anynull,
	  &status);
    indx++;
  }


  /* Get the correct value for NDumps by finding the number of HDUs and calculting
     it according to the version number of the ASP FITS file */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  if(!strcmp(hdr->gen.HdrVer,"Ver1.0")){
    hdr->redn.RNTimeDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				   many non-data tables we will be using */
  } 
  else if(!strcmp(hdr->gen.HdrVer,"Ver1.0.1")){
    hdr->redn.RNTimeDumps = (NumHDU-3)/2;
  }
  else{
    printf("Do not recognize ASP FITS file version number in header.\n");
    printf("This header %s. Exiting...\n",hdr->gen.HdrVer);fflush(stdout);
    exit(3);
  }
  /* Now check last dump to make sure it wrote properly if scan 
     was ctrl-c'd -- if so, reduce NDumps by 1 to avoid disaster */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  fits_movabs_hdu(Fin,NumHDU, &hdutype, &status);
  /* find NPtsProf */
  fits_get_num_rows(Fin, &NPtsProf, &status);status=0;    
  /* Won't happen for PSRFITS data this way */
  if ((hdr->redn.RNBinTimeDump!=(int)NPtsProf)){
    printf("Warning:  last dump did not write cleanly.\n");
    printf("Reducing number of dumps read in file by one...\n\n");
    fflush(stdout);
    hdr->redn.RNTimeDumps--;
    //RunMode.NDumps--; 
  }
  
  /* Fill in TDump value by finding the time stamp difference between the first 
     and second scans */
  if (hdr->redn.RNTimeDumps > 1) {
    if(!strcasecmp(hdr->gen.HdrVer,"Ver1.0")) {
      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT0", 0, &status);
      fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", 
		    &(DumpMiddleSecs), NULL, &status); status = 0;
    }
    else if (!strcasecmp(hdr->gen.HdrVer,"Ver1.0.1")) {
      fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF0", 0, &status);
      fits_read_key(Fin, TDOUBLE, "MIDSECS", &(DumpMiddleSecs), 
		    NULL, &status); status = 0;
    }
    hdr->redn.TDump = DumpMiddleSecs;
    /* Move to second data table , get central time stamp */
    if(!strcasecmp(hdr->gen.HdrVer,"Ver1.0")) {
      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT1", 0, &status);
      fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", 
		    &(DumpMiddleSecs),NULL, &status); status = 0;
    }
    else if (!strcasecmp(hdr->gen.HdrVer,"Ver1.0.1")) {
      fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF1", 0, &status);
      fits_read_key(Fin, TDOUBLE, "MIDSECS", &(DumpMiddleSecs), 
		    NULL, &status); status = 0;
    }
    if(DumpMiddleSecs < hdr->redn.TDump) DumpMiddleSecs += 86400.;
    
    hdr->redn.TDump = DumpMiddleSecs - hdr->redn.TDump;
  }
  else {
    /* We may have switched over to the next MJD, so take care of 
       this if it's the case */
    if (DumpMiddleSecs < (double)hdr->obs.StartTime) DumpMiddleSecs += 86400.;
    
    /* TDump is twice the difference between the middle time stamp 
       and the start time */
    hdr->redn.TDump = 2.* ((double)floor(DumpMiddleSecs) - (double)hdr->obs.StartTime);
  }
  

  /*  for (i = 0; i < 10; i++)
    {
      free(ttype[i]);
      free(tform[i]);
      free(tunit[i]);
      }  */

  retval = 0;

  return retval;
}
