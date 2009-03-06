/* ========================================================================= */
/* To write the ASP header, given the structure ASPHdr. It writes the header */
/* in FITS format                                                            */
/*                                                                           */
/* INPUTS : ASPHdr  : structure of header.                                   */
/*          Fin     : File Pointer, in "fitsfile" format.                    */
/*          FileName: char[]. Filename.                                      */
/*                                                                           */
/* R. Ramachandran, 15-March-2003, Berkeley.                                 */
/* ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPHeader.h"
#include "fitsio.h"

int WrtASPHdr(struct      ASPHdr *hdr,
	       fitsfile   *Fin)
{
  int     status=0, ival, i, ii, colnum, indx, retval=-1;
  char    *tunit[4*NCHMAX+10], *tform[4*NCHMAX+10], *ttype[4*NCHMAX+10], *junkch;

  junkch   = (char *) calloc(20, sizeof(char));


  ival = 1; fits_update_key(Fin, TLOGICAL, "SIMPLE", &ival, NULL, &status); status = 0;
  ival = 32; fits_update_key(Fin, TINT, "BITPIX", &ival, "FITS COMPATIBLE", &status); status = 0;
  ival = 0; fits_update_key(Fin, TINT, "NAXIS", &ival, "FITS COMPATIBLE", &status); status = 0;

  /* Scan Name of this scan */
  ffpky(Fin, TSTRING, "SCANNAME", hdr->gen.ScanName, NULL, &status); status = 0;  
  /* Ver of software used for analysis */
  ffpky(Fin, TSTRING, "SOFT_VER",  &hdr->gen.SoftVer,  NULL, &status); status = 0; 
  /* Computer platform used */  
  ffpky(Fin, TSTRING, "PLATFORM", &hdr->gen.Platform, NULL, &status); status = 0;
  /* Comments by Operator/Observer */
  ffpky(Fin, TSTRING, "COMMOPER", &hdr->gen.CommentOper, NULL, &status);status = 0;
  /* Header version no. */
  ffpky(Fin, TSTRING, "HDRVER", &hdr->gen.HdrVer, NULL, &status); status = 0; 
  /* Observer name */
  ffpky(Fin, TSTRING, "OBSERVER", &hdr->gen.Observer, NULL, &status); status = 0; 
  /* Project ID */
  ffpky(Fin, TSTRING, "PROJID", &hdr->gen.ProjID, NULL, &status); status = 0; 
  /* Observatory code */
  ffpky(Fin, TSTRING, "OBSVTY", &hdr->obs.ObsvtyCode, NULL, &status); status = 0; 
  /* Frontend name GREG/CH etc */
  ffpky(Fin, TSTRING, "FRONTEND", &hdr->gen.FEName, NULL, &status); status = 0; 
  /* Polarisation type CIRC/LIN */
  ffpky(Fin, TSTRING, "FD_POLN", &hdr->gen.FEPol, NULL, &status); status = 0; 
  /* Backend name ASP/ABPP/WAPP etc */
  ffpky(Fin, TSTRING, "BACKEND", &hdr->gen.BEName, NULL, &status); status = 0; 
  /* Backend Config filename */
  ffpky(Fin, TSTRING, "BECONFIG", &hdr->gen.BEConfFile, NULL, &status); status = 0;
  /* Observation mode PSR/CAL/SPECTRAL etc */
  ffpky(Fin, TSTRING, "OBS_MODE", &hdr->gen.ObsMode, NULL, &status); status = 0; 
  /* Source Name */
  ffpky(Fin, TSTRING, "SRC_NAME", &hdr->target.PSRName, NULL, &status); status = 0;
  /* Right Assention hh.hhhhhhhhhhhhhh */
  ffpky(Fin, TDOUBLE, "RA",      &hdr->target.RA, NULL, &status); status = 0; 
  /* Declination dd.dddddddddddddd */
  ffpky(Fin, TDOUBLE, "DEC",     &hdr->target.Dec, NULL, &status); status = 0; 
  /* Epoch of the coordinates */
  ffpky(Fin, TFLOAT, "EPOCH",   &hdr->target.Epoch, NULL, &status); status = 0;
  /* Coordinate mode for following 4 entries: J2000/GAL/ECLIPTIC/AZEL etc. */
  ffpky(Fin, TSTRING, "COORD_MD", &hdr->target.CoordMode, NULL,&status);status = 0;
  /* coordinate 1 at start time */
  ffpky(Fin, TDOUBLE, "STT_CRD1",  &hdr->target.StartCrd1,NULL,&status);status = 0;
  /* coordinate 2 at start time */
  ffpky(Fin, TDOUBLE, "STT_CRD2",&hdr->target.StartCrd2, NULL, &status);status = 0;
  /* coordinate 1 at stop time */
  ffpky(Fin, TDOUBLE, "STP_CRD1",&hdr->target.StopCrd1, NULL, &status); status = 0;
  /* coordinate 2 at stop time */
  ffpky(Fin, TDOUBLE, "STP_CRD2",&hdr->target.StopCrd2, NULL, &status);status = 0; 
  /* Track mode TRACK/SLEW/DRIFT etc. */
  ffpky(Fin, TSTRING, "TRK_MODE",&hdr->target.TrackMode, NULL, &status);status = 0;
  /* Observation length (s) */
  ffpky(Fin, TDOUBLE, "SCANLEN", &hdr->obs.ObsLength, NULL, &status); status = 0; 
  /* Start Date/Month/Year */
  ffpky(Fin, TSTRING, "STT_DATE", &hdr->obs.StartDate, NULL, &status); status = 0; 
  /* Start Time in UT */
  ffpky(Fin, TSTRING, "STT_TIME", &hdr->obs.StartUT, NULL, &status); status = 0; 
  /* Integer MJD of starting sample */
  ffpky(Fin, TINT,    "STT_IMJD", &hdr->obs.IMJDStart, NULL, &status); status = 0; 
  /* Start Time (INT) past 0.0h */
  ffpky(Fin, TINT, "STT_SMJD", &hdr->obs.StartTime, NULL, &status); status = 0; 
  /* Fractional second of start time */
  ffpky(Fin, TDOUBLE, "STT_FRAC", &hdr->obs.StartFSec, NULL, &status); status = 0; 
  /* Clock offset (sec) */
  ffpky(Fin, TDOUBLE, "STT_OFFS", &hdr->obs.ClockOffset, NULL, &status);status = 0;
  /* Start LST */
  ffpky(Fin, TDOUBLE, "STT_LST", &hdr->target.StartLST, NULL, &status); status = 0;
  /* Start LST */
  ffpky(Fin, TFLOAT, "ION_RM", &hdr->obs.IonRM, NULL, &status); status = 0; 
  /* Sampling interval (sec) */
  ffpky(Fin, TDOUBLE, "SAMP_INT", &hdr->obs.SampInterval, NULL,&status);status = 0;
  /* Output data type String (IQUVLT etc) */
  ffpky(Fin, TSTRING, "OP_STRNG", hdr->obs.OPString, NULL, &status); status = 0; 
  /* Out put scale 1-Arbit; 2-Jy; 3-K;  */
  ffpky(Fin, TINT, "OP_SCALE", &hdr->obs.OPScale, NULL, &status); status = 0; 
  /* Dynamic DC subtracted? 0/1 */
  ffpky(Fin, TINT, "DYN_DC", &hdr->obs.DynDC, NULL, &status); status = 0; 
  /* Dynamic DC correction time scale (s) */
  ffpky(Fin, TFLOAT, "DYN_TIME", &hdr->obs.DynTime, NULL, &status); status = 0; 
  /* No. of bits per sample */
  ffpky(Fin, TINT, "NBITS", &hdr->obs.NBitPerSamp, NULL, &status); status = 0; 
  /* Interference corrected? 0/1 */
  ffpky(Fin, TINT, "INTERF", &hdr->obs.InterfString, NULL, &status); status = 0; 
  /* Zap filename */
  ffpky(Fin, TSTRING, "ZAPFILE", &hdr->obs.ZapFileName, NULL, &status); status = 0;
  /* Number of bins per profile */
  ffpky(Fin, TINT, "NPTSPROF", &hdr->redn.RNBinTimeDump, NULL,&status); status = 0;
  /* Number of bins per profile */
  ffpky(Fin, TINT, "NDUMPS", &hdr->redn.RNTimeDumps, NULL, &status); status = 0; 
  /* Overall Central Sky Frequency */
  ffpky(Fin, TDOUBLE, "FSKYCENT", &hdr->obs.FSkyCent, NULL, &status); status = 0; 
  /* Forward/Reverse band? [1 -> forward / -1 -> reverse] */
  ffpky(Fin, TINT, "SIDEBAND", &hdr->obs.Sideband, NULL, &status); status = 0; 
  
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
    strcpy(&junkch[0], "CFRQ\0"); sprintf(&junkch[4], "%d", i);
    strcpy(ttype[++colnum], junkch);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch[0], "BW\0"); sprintf(&junkch[2], "%d", i);
    strcpy(ttype[++colnum], junkch);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch[0], "TSYS\0"); sprintf(&junkch[4], "%d", i);
    strcpy(ttype[++colnum], junkch);
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch[0], "TCAL\0"); sprintf(&junkch[4], "%d", i);
    strcpy(ttype[++colnum], junkch);
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

  fits_create_tbl(Fin, ASCII_TBL, 0, colnum, ttype, tform, tunit, "BECONFIG", &status);

  indx = 1;
  fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->obs.NChan, &status); indx++;

  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TDOUBLE, indx, 1, 1, 1, &hdr->obs.ChanFreq[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TDOUBLE, indx, 1, 1, 1, &hdr->obs.ChanWidth[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TFLOAT, indx, 1, 1, 1, &hdr->obs.ChanTSys[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TFLOAT, indx, 1, 1, 1, &hdr->obs.ChanTCal[i], &status);
    indx++;
  }
  fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->asp.SplitFac, &status);
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
    strcpy(&junkch[0], "NCHIRP\0"); sprintf(&junkch[6], "%d", i); strcpy(ttype[++colnum], junkch); }
  for (i=0; i<hdr->obs.NChan; i++) {
    strcpy(&junkch[0], "NOVRLAP\0"); sprintf(&junkch[7], "%d", i); strcpy(ttype[++colnum], junkch); }

  colnum  = -1;
  strcpy(tunit[++colnum], "pc/cc");
  strcpy(tunit[++colnum], "");
  strcpy(tunit[++colnum], "rad/m^2");
  strcpy(tunit[++colnum], "");
  for (i=0; i<hdr->obs.NChan; i++)
    strcpy(tunit[++colnum], "");
  for (i=0; i<hdr->obs.NChan; i++)
    strcpy(tunit[++colnum], "");


  fits_create_tbl(Fin, ASCII_TBL, 0, ++colnum, ttype, tform, tunit, "COHDDISP", &status);
  //printf("fits_create_tbl status = %d\n", status);

  indx = 1;
  fits_write_col(Fin, TDOUBLE, indx, 1, 1, 1, &hdr->obs.DM, &status); indx++;
  fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->obs.DMMethod, &status); indx++;
  fits_write_col(Fin, TFLOAT, indx, 1, 1, 1, &hdr->obs.RM, &status); indx++;
  fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->obs.RMMethod, &status); indx++;
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->obs.ChanChirpLen[i], &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    fits_write_col(Fin, TINT, indx, 1, 1, 1, &hdr->obs.ChanOverlap[i], &status);
    indx++;
  }




  retval = 0;
  return retval;

}
