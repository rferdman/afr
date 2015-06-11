/* ========================================================================= */
/* Write output PSRFITS header, and set up tables for output data.           */
/*                                                                           */
/* INPUTS : ASPHdr  : structure of header.                                   */
/*          Fin     : File Pointer, in "fitsfile" format.                    */
/*                                                                           */
/* R. Ferdman, 4-March-2013, McGill University, Montreal.                    */
/* ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "ASPHeader.h"
#include "fitsio.h"

int WrtPSRFITSHdr(struct    ASPHdr *hdr,
	          fitsfile  *Fout)
{
  int     status=0, ival, i, ii, colnum, indx, retval=-1;
  double  ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec;
  // char    *tunit[4*NCHMAX+10], *tform[4*NCHMAX+10], *ttype[4*NCHMAX+10];
  char    junkch1[2048], junkch2[2048], junkch3[2048];
  time_t  cur_time;
  struct tm *utctime;

  //  junkch1   = (char *) calloc(20, sizeof(char));


  /*     ival = 1; fits_update_key(Fout, TLOGICAL, "SIMPLE", &ival, NULL, &status); status = 0;
	 ival = 32; fits_update_key(Fout, TINT, "BITPIX", &ival, "FITS COMPATIBLE", &status); status = 0;
	 ival = 0; fits_update_key(Fout, TINT, "NAXIS", &ival, "FITS COMPATIBLE", &status); status = 0;     */

  /* Fill in primary header with current header variables, and change what 
     is needed */

  /* Create primary header */
  
  if(fits_create_img(Fout, 8, 0, NULL, &status)  != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot create primary table.\n");
    return -1;
  }


  /* Start adding keywords */

  /* Header version */
  sprintf(junkch1, "5.0");
  if(fits_write_key(Fout, TSTRING, "HDRVER", junkch1, "Header version", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword HDRVER.\n");
    return -1;
  }
 
  /* Type of FITS file */
  sprintf(dummy_txt, "PSRFITS");
  if(fits_write_key(Fout, TSTRING, "FITSTYPE", junkch1, "FITS definition for pulsar data files", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FITSTYPE.\n");
    return -1;
  }
  
  /* File creation date (YYYY-MM-DDThh:mm:ss UTC) */ 
  /* Get current time */
  curtime = time(NULL);
  /* Convert time to gmtine representation */
  utctime = gmtime(&curtime);
  /* String it up */
  strftime(junkch1, 256, "%Y-%m-%dT%H:%M:%S");
  //  sprintf(junkch1, "%sT%s", hdr->filecont.FileDate, hdr->filecont.FileUT);
  if(fits_write_key(Fout, TSTRING, "DATE", junkch1, "File creation date (YYYY-MM-DDThh:mm:ss UTC)", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword DATE.\n");
    return -1;
  }
  
  /* Observer name */
  if(fits_write_key(Fout, TSTRING, "OBSERVER", hdr->gen.Observer, "Observer name(s)", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword OBSERVER.\n");
    return -1;
  }

/* Project ID */
  if(fits_write_key(Fout, TSTRING, "PROJID", hdr->gen.ProjID, "Project name", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword PROJID.\n");
    return -1;
  } 

  /* Telescope identifier */
  if(fits_write_key(Fout, TSTRING, "TELESCOP", hdr->obs.ObsvtyCode, "Telescope name", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword TELESCOP.\n");
    return -1;
  }
 
/* Antenna ITRF X, Y, and Z coordinates */
  if(fits_write_key(Fout, TDOUBLE, "ANT_X", &hdr->target.Ant_X, "Antenna ITRF X-coordinate", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword ANT_X.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "ANT_X", "m", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword ANT_X.\n");
    return -1;
  }

  if(fits_write_key(Fout, TDOUBLE, "ANT_Y", &hdr->target.Ant_Y, "Antenna ITRF Y-coordinate", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword ANT_Y.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "ANT_Y", "m", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword ANT_Y.\n");
    return -1;
  }

  if(fits_write_key(Fout, TDOUBLE, "ANT_Z", &hdr->target.Ant_Z, "Antenna ITRF Z-coordinate", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword ANT_Z.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "ANT_Z", "m", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword ANT_Z.\n");
    return -1;
  }


/* Frontend name GREG/CH etc */
  if(fits_write_key(Fout, TSTRING, "FRONTEND", hdr->gen.FEName, "Receiver ID", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FRONTEND.\n");
    return -1;
  }

/* Beam ID for multiple-beam systems */
  if(fits_write_key(Fout, TSTRING, "IBEAM", hdr->gen.IBeam, "Beam ID for multibeam systems", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword IBEAM.\n");
    return -1;
  }
  
/* Number of receiver polarisation channels */
  if(fits_write_key(Fout, TSTRING, "NRCVR", hdr->gen.NRcvr, "Number of receiver polarisation channels", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword NRCVR.\n");    
    return -1;
  } 

/* Polarisation type CIRC/LIN */
  if (!strncmp(hdr->gen.FEPol, "L", 1))
    strcpy(junkch1, "LIN");
  else if (!strncmp(hdr->gen.FEPol, "C", 1))
    strcpy(junkch1, "CIRC");
  else {
    fprintf(stderr,"WrtPSRFITSHdr ERROR: Unrecognized polarization basis %s\n",
	    hdr->gen.FEPol);
    return -1;
  }
  if(fits_write_key(Fout, TSTRING, "FD_POLN", hdr->gen.FEPol, "Polarisation type CIRC/LIN", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_POLN.\n");    
    return -1;  
  } 


/* A/B polarization swap -- +/- 1.  +1 is LIN: A=X, B=Y, CIRC:A=L, B=R */
  if(fits_write_key(Fout, TINT, "FD_HAND", &hdr->gen.FDHand, "+/- 1.  +1 is LIN: A=X, B=Y, CIRC:A=L, B=R", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_HAND.\n");    
    return -1;
  }
/* FA of E vector for equal sig in A&B */
  if(fits_write_key(Fout, TFLOAT, "FD_SANG", &hdr->gen.FDSang, "FA of E vector for equal sig in A&B", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_SANG.\n");   
    return -1; 
  }
  if(fits_write_key_unit(Fout, "FD_SANG", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword FD_SANG.\n");
    return -1;
  }

/* Phase of A*B for injected cal */ 
  if(fits_write_key(Fout, TFLOAT, "FD_XYPH", &hdr->gen.FDXYph, "Phase of A^* B for injected cal", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_XYPH.\n");   
    return -1;
  }
  if(fits_write_key_unit(Fout, "FD_XYPH", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword FD_XYPH.\n");
    return -1;
  }


  /* Name of backend instrument used */
  if(fits_write_key(Fout, TSTRING, "BACKEND", hdr->gen.BEName, "Backend ID", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BACKEND.\n");
    return -1;
  }

/* Backend Config filename */
  if(fits_write_key(Fout, TSTRING, "BECONFIG", hdr->gen.BEConfFile, "Backend configuration file name", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BECONFIG.\n");
    return -1;
  }
  
/* 0/+1/-1 BE cross-phase: 0 unknown, +/-1 std/rev */
  if(fits_write_key(Fout, TINT, "BE_PHASE", &hdr->gen.BEPhase, "0/+1/-1 BE cross-phase: 0 unknown, +/-1 std/rev", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BE_PHASE.\n");
    return -1;
  }

/* 0/1 BE downconversion conjugation */
  if(fits_write_key(Fout, TINT, "BE_DCC", &hdr->gen.BEDownConv, "0/1 BE downconversion conjugation", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BE_DCC.\n");
    return -1;
  } 

/* Backend propn delay from digitiser in seconds */
  if(fits_write_key(Fout, TFLOAT, "BE_DELAY", &hdr->gen.BEDelay, "Backend propn delay from digitiser in seconds", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BE_DELAY.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "BE_DELAY", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword BE_DELAY.\n");
    return -1;
  }

/* On-line cycle time */
  if(fits_write_key(Fout, TDOUBLE, "TCYCLE", &hdr->gen.TCycle, "On-line cycle time", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword TCYCLE.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "TCYCLE", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword TCYCLE.\n");
    return -1;
  }

  /* Mode of observing (usually PSR or CAL) */
  if(fits_write_key(Fout, TSTRING, "OBS_MODE", hdr->gen.ObsMode, "", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword OBS_MODE.\n");
    return -1;
  }

/* Date of observation (YYYY-MM-DDThh:mm:ss UTC) */ 
  sprintf(junkch1, "%sT%s", hdr->obs.StartDate, hdr->obs.StartUT);
  if(fits_write_key(Fout, TSTRING, "DATE-OBS", tempstr1, "Date of observation (YYYY-MM-DDThh:mm:ss UTC)", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword DATE-OBS.\n");
    return -1;
  }

   /* Centre observing frequency */
  if(fits_write_key(Fout, TDOUBLE, "OBSFREQ", hdr->obs.FSkyCent, "Centre frequency for observation", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword OBSFREQ.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "OBSFREQ", "MHz", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword OBSFREQ.\n");
    return -1;
  }

  /* Bandwidth for observation [MHz] */
  /* Multiply by Sideband to get it in proper format for psrfits definition */
  hdr->obs.BW *= (double)hdr->obs.Sideband;
  if(fits_write_key(Fout, TDOUBLE, "OBSBW", &hdr->obs.BW, "Bandwidth for observation", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword OBSBW.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "OBSBW", "MHz", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword OBSBW.\n");
    return -1;
  }

/* Number of frequency channels (original number) */
  if(fits_write_key(Fout, TINT, "OBSNCHAN", &hdr->obs.NChanOrig, "Number of frequency channels (original number)", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword OBSNCHAN.\n");
    return -1;
  }


  /* DM used for online dedispaersion [cm-3 pc] */
  if(fits_write_key(Fout, TDOUBLE, "CHAN_DM", &hdr->obs.ChanDM, "DM used for online dedispaersion", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CHAN_DM.\n");
    return -1;
  } 
  if(fits_write_key_unit(Fout, "CHAN_DM", "cm-3 pc", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword CHAN_DM.\n");
    return -1;
  }
 
  
  /* Name or ID for pointing centre (multibeam feeds) */
  if(fits_write_key(Fout, TSTRING, "PNT_ID", hdr->gen.PntID, "Name or ID for pointing centre (multibeam feeds)", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword PNT_ID.\n");
    return -1;
  }
  
  
  /* Pulsar or source name */
  if(fits_write_key(Fout, TSTRING, "SRC_NAME", hdr->target.PSRName, "Source or scan ID", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword SRC_NAME.\n");
    return -1;
  }
  
  /*===========================================*/

/* Coordinate mode for following 4 entries  J2000/GAL/ECLIPTIC/AZEL etc. */
if(fits_write_key(Fout, TSTRING, "COORD_MD", hdr->target.CoordMode, "Coordinate mode (J2000/GAL/ECLIPTIC/AZEL etc.)", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword COORD_MD.\n");
    return -1;
 }

/* Epoch of the coordinates */
/* Read in as string and convert depending on PSRFITS version */
  if(fits_write_key(Fout, TFLOAT, "EQUINOX", &hdr->target.Epoch, "Epoch of the coordinates", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword EQUINOX.\n");
    return -1; 
  }
 
/* Right Assention hh:mm:ss.sss */
  ra_hr = floor(hdr->target.RA);
  ra_min = floor((hdr->target.RA - ra_hr)*60.);
  ra_sec = ((hdr->target.RA - ra_hr)*60. - ra_min) * 60.;
  sprintf(junkch1, "%2d:%2d:%lf", (int)ra_hr, (int)ra_min, ra_sec);
  if(fits_write_key(Fout, TSTRING, "RA", junkch1, "Right Assention (hh:mm:ss.sss)", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword RA.\n");
    return -1;
  }

/* Declination -dd:mm:ss.sss */
  dec_deg = floor(hdr->target.Dec);
  dec_min = floor((hdr->target.Dec - dec_deg)*60.);
  dec_sec = ((hdr->target.Dec - dec_deg)*60. - dec_min) * 60.;
  if(hdr->target.Dec < 0): 
    sprintf(junkch1, "-%2d:%2d:%lf", (int)dec_deg, (int)dec_min, dec_sec);
  else
    sprintf(junkch1, "+%2d:%2d:%lf", (int)dec_deg, (int)dec_min, dec_sec);
  if(fits_write_key(Fout, TSTRING, "DEC", junkch1, "Declination (-dd:mm:ss.sss)", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword DEC.\n");
    return -1;
  }


  if(deg > 0.) sgn = 1.; else sgn = -1;
  hdr->target.Dec = deg + sgn*(min/60. + sec/3600.);

/* Beam major axis length */
  if(fits_write_key(Fout, TDOUBLE, "BMAJ", &hdr->target.BeamMajor, "Beam major axis length", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BMAJ.\n");
    return -1;  
  }
  if(fits_write_key_unit(Fout, "BMAJ", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword BMAJ.\n");
    return -1;
  }
 
/* Beam minor axis length */
  if(fits_write_key(Fout, TDOUBLE, "BMIN", &hdr->target.BeamMinor, "Beam minor axis length", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BMIN.\n");
    return -1;   
  }
  if(fits_write_key_unit(Fout, "BMIN", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword BMIN.\n");
    return -1;
  }

/* Beam position angle */
  if(fits_write_key(Fout, TDOUBLE, "BPA", &hdr->target.BeamPA, "Beam position angle", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword BPA.\n");
    return -1;   
  }
  if(fits_write_key_unit(Fout, "BPA", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword BPA.\n");
    return -1;
  }

/* coordinate 1 at start time */
  if(fits_write_key(Fout, TSTRING, "STT_CRD1",hdr->target.StartCrd1, "Coordinate 1 at start time", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_CRD1.\n");
    return -1; 
  }

/* coordinate 2 at start time */
  if(fits_write_key(Fout, TSTRING, "STT_CRD2",hdr->target.StartCrd2, "Coordinate 2 at start tim", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_CRD2.\n");
    return -1; 
  }

/* coordinate 1 at stop time */
  if(fits_write_key(Fout, TSTRING, "STP_CRD1", hdr->target.StopCrd1, "Coordinate 1 at stop time", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STP_CRD1.\n");
    return -1; 
  }

/* coordinate 2 at stop time */
  if(fits_write_key(Fout, TSTRING, "STP_CRD2", hdr->target.StopCrd2, "Coordinate 2 at stop time", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STP_CRD2.\n");
    return -1; 
  }
  
/* Track mode TRACK/SLEW/DRIFT/SCANGO/SCANLAT etc. */
  if(fits_write_key(Fout, TSTRING, "TRK_MODE", hdr->target.TrackMode, "Track mode (TRACK/SLEW/DRIFT/SCANGO/SCANLAT etc.)", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword TRK_MODE.\n");
    return -1; 
  }

/* Requested observation length (s) */
  if(fits_write_key(Fout, TDOUBLE, "SCANLEN", &hdr->obs.ObsLength, "Requested observation length", &status) !=0 ){ 
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword SCANLEN.\n");
    return -1; 
  }
  if(fits_write_key_unit(Fout, "SCANLEN", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword SCANLEN.\n");
    return -1;
  }

/* Feed track mode -- FA/CPA/SPA/TPA */
  if(fits_write_key(Fout, TSTRING, "FD_MODE", hdr->gen.FeedMode, "Feed track mode -- FA/CPA/SPA/TPA", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_MODE.\n");
    return -1; 
  }

/* Feed/Position angle requested */
  if(fits_write_key(Fout, TDOUBLE, "FA_REQ", &hdr->gen.FAReq, "Feed/Position angle requested", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword FD_REQ.\n");
    return -1; 
  }
  if(fits_write_key_unit(Fout, "FA_REQ", "deg", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword FA_REQ.\n");
    return -1;
  }

/* Cal mode (OFF, SYNC, EXT1, EXT2) */
  if(fits_write_key(Fout, TSTRING, "CAL_MODE", hdr->gen.CalMode, "Cal mode (OFF, SYNC, EXT1, EXT2)", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CAL_MODE.\n");
    return -1; 
  }

/* Cal modulation frequency in Hz */  
  if(fits_write_key(Fout, TDOUBLE, "CAL_FREQ", &hdr->obs.CalFreq, "Cal modulation frequency", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CAL_FREQ.\n");
    return -1; 
  }
  if(fits_write_key_unit(Fout, "CAL_FREQ", "Hz", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword CAL_FREQ.\n");
    return -1;
  }


/* Cal duty cycle */  
  if(fits_write_key(Fout, TDOUBLE, "CAL_DCYCLE", &hdr->obs.CalDutyCycle, "Cal duty cycle", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CAL_DCYCLE.\n");
    return -1; 
  }

/* Cal phase (wrt start time) */
  if(fits_write_key(Fout, TDOUBLE, "CAL_PHS", &hdr->obs.CalPhase, "Cal phase (wrt start time)", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CAL_PHS.\n");
    return -1; 
  }

/* Number of states in cal pulse */
  if(fits_write_key(Fout, TDOUBLE, "CAL_NPHS", &hdr->obs.CalNStates, "Number of states in cal pulse", &status) !=0 ){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword CAL_NPHS.\n");
    return -1; 
  }    
  
  /*==========================================*/

  /* Integer part of start MJD */
  if(fits_write_key(Fout, TINT, "STT_IMJD", hdr->obs.IMJDStart, "Integer part of start MJD", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_IMJD.\n");
    return -1;
  }

  /* Start time in number of seconds since 00:00 UT */
  if(fits_write_key(Fout, TINT, "STT_SMJD", hdr->obs.StartTime, "Start time in number of seconds past 00:00 UT", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_SMJD.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "STT_SMJD", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword STT_SMJD.\n");
    return -1;
  }


  /* Start time of offset (sec) */
  if(fits_write_key(Fout, TDOUBLE, "STT_OFFS", &hdr->obs.ClockOffset, "Start time of offset", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_OFFS.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "STT_OFFS", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword STT_OFFS.\n");
    return -1;
  }

  /* Start LST */
  if(fits_write_key(Fout, TDOUBLE, "STT_LST", &hdr->target.StartLST, "Start LST", &status) != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write keyword STT_LST.\n");
    return -1;
  }
  if(fits_write_key_unit(Fout, "STT_LST", "s", &status) != 0){
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot write units for keyword STT_LST.\n");
    return -1;
  }



  /***********  DONE PRIMARY HEADER  **************/

  /***********  WRITE OTHER TABLES WE MAY USE LATER ON ************/

  /* Will replace NULLs, etc. here with useful stuff in the future */
  if(fits_create_tbl(Fout, BINARY_TBL, 0, 0, NULL, NULL, NULL, "HISTORY", &status)  != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot create table HISTORY.\n");
    return -1;
  }
  
  /* History table will include the AFR command line used to create the output data file, 
     plus other relevant info */

  


  /***********  CREATE SUBINT HEADER AND SET UP TABLE FOR WRITING  ************/

  
  /* Set up column names, format, and units */
  char *ttype[] = {"INDEXVAL", "TSUBINT", "OFFS_SUB", "LST_SUB", "RA_SUB", "DEC_SUB", "GLON_SUB", "GLAT_SUB", "FD_ANG", "POS_ANG", "PAR_ANG", "TEL_AZ", "TEL_ZEN", "AUX_DM", "AUX_RM", "DAT_FREQ", "DAT_WTS", "DAT_OFFS", "DAT_SCL", "DATA"};
  sprintf(junkch1, "%ldE", hdr->obs.NChan);
  sprintf(junkch2, "%ldE", hdr->obs.NChan*hdr->obs.NPoln);
  sprintf(junkch3, "%ldI", hdr->obs.NChan*hdr->obs.NPoln*hdr->redn.RNBinTimeDump)
  char *tform[] = {"1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1E", "1E", "1E", "1E", "1E", "1D", "1D", junkch1, junkch1, junkch2, junkch2, junkch3};
  char *tunit[] = {"", "s", "s", "s", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "cm-3 pc", "rad m-2", "MHz", "", "", "", ""};

  /* Create the table */
  /* Hardcoding 20 columns for now */
  if(fits_create_tbl(Fin, BINARY_TBL, hdr->redn.RNTimeDumps, 20, ttype, tform, tunit, "SUBINT", &status)  != 0) {
    fprintf(stderr, "WrtPSRFITSHdr ERROR: Cannot create table SUBINT.\n");
    return -1;
  }

  /* Start adding keywords */
  






  /* ======================================================================= */
  //  printf("Going to create ASCII table\n");

  if (hdr->obs.NChan == 0) {
    printf("\nNumber of ASP channels == zero!!!\n");
    printf("   .... How is it POSSIBLE?!?!!\n");
    exit(0);
  }

  for (ii = 0; ii < (4*(hdr->obs.NChan)+10); ii++)
    {
      ttype[ii] = (char *) malloc(64);
      tform[ii] = (char *) malloc(64);
      tunit[ii] = (char *) malloc(64);
    }

  //printf("Arrays initialised  %d\n", (4* hdr->obs.NChan + 10));

  colnum = -1;
  strcpy(tform[++colnum], "I3");
  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "F10.5"); }
  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "F10.5"); }
  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "F7.1"); }
  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "F7.2"); }
  strcpy(tform[++colnum], "I3");
  //printf("Units defined\n");

  colnum  = -1;
  strcpy(ttype[++colnum], "NChan");
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "CFRQ\0"); sprintf(&junkch1[4], "%d", i);
    strcpy(ttype[++colnum], junkch1);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "BW\0"); sprintf(&junkch1[2], "%d", i);
    strcpy(ttype[++colnum], junkch1);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "TSYS\0"); sprintf(&junkch1[4], "%d", i);
    strcpy(ttype[++colnum], junkch1);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "TCAL\0"); sprintf(&junkch1[4], "%d", i);
    strcpy(ttype[++colnum], junkch1);
  }
  strcpy(ttype[++colnum], "SPLITFAC");
  //printf("Labels defined\n");

  colnum = -1;
  strcpy(tunit[++colnum], "");
  for (i=0; i<hdr->obs.NChan; i++) strcpy(tunit[++colnum], "MHz");
  for (i=0; i<hdr->obs.NChan; i++) strcpy(tunit[++colnum], "MHz");
  for (i=0; i<hdr->obs.NChan; i++) strcpy(tunit[++colnum], "K");
  for (i=0; i<hdr->obs.NChan; i++) strcpy(tunit[++colnum], "K");
  strcpy(tunit[++colnum], "");

  colnum++;
  //printf("COLNUM = %d\n", colnum);

  fits_create_tbl(Fout, ASCII_TBL, 0, colnum, ttype, tform, tunit, "BECONFIG", &status);

  indx = 1;
  fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->obs.NChan, &status); indx++;

  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TDOUBLE, indx, 1, 1, 1, &hdr->obs.ChanFreq[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TDOUBLE, indx, 1, 1, 1, &hdr->obs.ChanWidth, &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TFLOAT, indx, 1, 1, 1, &hdr->obs.ChanTSys[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TFLOAT, indx, 1, 1, 1, &hdr->obs.ChanTCal[i], &status);
    indx++;
  }
  fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->asp.SplitFac, &status);
  //printf("SPLITFAC = %d\n",hdr->asp.SplitFac);fflush(stdout);

  /* ======================================================================= */
  /* Dedispersion and Rotation Measure table */

  colnum = -1;
  strcpy(tform[++colnum], "F10.5");    /* Dispersion Measure (pc/cc) */
  strcpy(tform[++colnum], "I1");       /* Dedispersion type 1 --> Coherent; 
					2 --> Incoherent; 3 --> Coh/Incoh */

  strcpy(tform[++colnum], "F10.5");    /* Rotation Measure (pc/cc) */
  strcpy(tform[++colnum], "I1");       /* Faraday Correction type 1 --> Coherent;
					2 --> Incoherent; 3 --> Coh/Incoh */

  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "I8"); }  /* Number of pts in each chirp function */
  for (i=0; i<hdr->obs.NChan; i++) { strcpy(tform[++colnum], "I8"); }  /* Number of overlapping points */

  colnum  = -1;
  strcpy(ttype[++colnum], "DM");
  strcpy(ttype[++colnum], "DDPTYPE");
  strcpy(ttype[++colnum], "RM");
  strcpy(ttype[++colnum], "DRMTYPE");  
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "NCHIRP\0"); sprintf(&junkch1[6], "%d", i); strcpy(ttype[++colnum], junkch1); }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch1[0], "NOVRLAP\0"); sprintf(&junkch1[7], "%d", i); strcpy(ttype[++colnum], junkch1); }

  colnum  = -1;
  strcpy(tunit[++colnum], "pc/cc");
  strcpy(tunit[++colnum], "");
  strcpy(tunit[++colnum], "rad/m^2");
  strcpy(tunit[++colnum], "");
  for (i=0; i<hdr->obs.NChan; i++)
    strcpy(tunit[++colnum], "");
  for (i=0; i<hdr->obs.NChan; i++)
    strcpy(tunit[++colnum], "");


  fits_create_tbl(Fout, ASCII_TBL, 0, ++colnum, ttype, tform, tunit, "COHDDISP", &status);
  //printf("fits_create_tbl status = %d\n", status);

  indx = 1;
  fits_write_col(Fout, TDOUBLE, indx, 1, 1, 1, &hdr->obs.DM, &status); indx++;
  fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->obs.DMMethod, &status); indx++;
  fits_write_col(Fout, TFLOAT, indx, 1, 1, 1, &hdr->obs.RM, &status); indx++;
  fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->obs.RMMethod, &status); indx++;
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->obs.ChanChirpLen[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fout, TINT, indx, 1, 1, 1, &hdr->obs.ChanOverlap[i], &status);
    indx++;
  }


  for (ii = 0; ii < (4*(hdr->obs.NChan)+10); ii++)
    {
      free(ttype[ii]);
      free(tform[ii]);
      free(tunit[ii]);
    }


  retval = 0;
  return retval;

}
