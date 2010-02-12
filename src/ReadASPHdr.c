#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPHeader.h"
#include "fitsio.h"

int ReadASPHdr(struct      ASPHdr *hdr,
		fitsfile   *Fin)
{
  int     i, status, knul=0, anynull, indx, retval;
  double  dnul=0.0;
  float   fnul=0;
  long    nrows=0;
  char    tblname[40], *ttype[10], *tform[10], *tunit[10];

  retval = -1;

  ffgky(Fin, TSTRING, "SCANNAME", hdr->gen.ScanName, NULL, &status); status=0;
/* Ver of software used for analysis */
  ffgky(Fin, TSTRING, "SOFT_VER", hdr->gen.SoftVer,  NULL, &status); status=0; 
/* Computer platform used */  
  ffgky(Fin, TSTRING, "PLATFORM", hdr->gen.Platform, NULL, &status); status=0;  
/* Comments by Operator/Observer */
  ffgky(Fin, TSTRING, "COMMOPER", hdr->gen.CommentOper, NULL, &status); status=0;  
/* Header version no. */
  ffgky(Fin, TSTRING, "HDRVER", hdr->gen.HdrVer, NULL, &status); status=0; 
/* Observer name */
  ffgky(Fin, TSTRING, "OBSERVER", hdr->gen.Observer, NULL, &status); status=0; 
/* Project ID */
  ffgky(Fin, TSTRING, "PROJID", hdr->gen.ProjID, NULL, &status); status=0; 
/* Observatory code */
  ffgky(Fin, TSTRING, "OBSVTY", hdr->obs.ObsvtyCode, NULL, &status); status=0; 
/* Frontend name GREG/CH etc */
  ffgky(Fin, TSTRING, "FRONTEND", hdr->gen.FEName, NULL, &status); status=0; 
/* Polarisation type CIRC/LIN */
  ffgky(Fin, TSTRING, "FD_POLN", hdr->gen.FEPol, NULL, &status); status=0; 
/* Backend name ASP/ABPP/WAPP etc */
  ffgky(Fin, TSTRING, "BACKEND", hdr->gen.BEName, NULL, &status); status=0; 
/* Backend Config filename */
  ffgky(Fin, TSTRING, "BECONFIG", hdr->gen.BEConfFile, NULL, &status); status=0; 
/* Observation mode PSR/CAL/SPECTRAL etc */
  ffgky(Fin, TSTRING, "OBS_MODE", hdr->gen.ObsMode, NULL, &status); status=0; 
/* Source Name */
  ffgky(Fin, TSTRING, "SRC_NAME", hdr->target.PSRName, NULL, &status); status=0; 
/* Right Assention hh.hhhhhhhhhhhhhh */
  ffgky(Fin, TDOUBLE, "RA", &hdr->target.RA, NULL, &status); status=0; 
/* Declination dd.dddddddddddddd */
  ffgky(Fin, TDOUBLE, "DEC", &hdr->target.Dec, NULL, &status); status=0; 
/* Epoch of the coordinates */
  ffgky(Fin, TFLOAT, "EPOCH", &hdr->target.Epoch, NULL, &status); status=0; 
/* Coordinate mode for following 4 entries  J2000/GAL/ECLIPTIC/AZEL etc. */
  ffgky(Fin, TSTRING, "COORD_MD", hdr->target.CoordMode, NULL, &status);status=0; 
/* coordinate 1 at start time */
  ffgky(Fin, TDOUBLE, "STT_CRD1",&hdr->target.StartCrd1, NULL, &status);status=0; 
/* coordinate 2 at start time */
  ffgky(Fin, TDOUBLE, "STT_CRD2",&hdr->target.StartCrd2, NULL, &status);status=0; 
/* coordinate 1 at stop time */
  ffgky(Fin, TDOUBLE, "STP_CRD1", &hdr->target.StopCrd1, NULL, &status);status=0; 
/* coordinate 2 at stop time */
  ffgky(Fin, TDOUBLE, "STP_CRD2", &hdr->target.StopCrd2, NULL, &status);status=0; 
/* Track mode TRACK/SLEW/DRIFT etc. */
  ffgky(Fin, TSTRING, "TRK_MODE", hdr->target.TrackMode, NULL, &status);status=0; 
/* Observation length (s) */
  ffgky(Fin, TDOUBLE, "SCANLEN", &hdr->obs.ObsLength, NULL, &status); status=0; 
/* Start Date/Month/Year */
  ffgky(Fin, TSTRING, "STT_DATE", hdr->obs.StartDate, NULL, &status); status=0; 
/* Start Time in UT */
  ffgky(Fin, TSTRING, "STT_TIME", hdr->obs.StartUT, NULL, &status); status=0; 
/* Integer MJD of starting sample */
  ffgky(Fin, TINT,    "STT_IMJD", &hdr->obs.IMJDStart, NULL, &status); status=0; 
/* Start Time (INT) past 0.0h */
  ffgky(Fin, TINT, "STT_SMJD", &hdr->obs.StartTime, NULL, &status); status=0; 
/* Fractional second of start time */
  ffgky(Fin, TDOUBLE, "STT_FRAC", &hdr->obs.StartFSec, NULL, &status); status=0; 
/* Clock offset (sec) */
  ffgky(Fin, TDOUBLE, "STT_OFFS", &hdr->obs.ClockOffset, NULL, &status);status=0; 
/* Start LST */
  ffgky(Fin, TDOUBLE, "STT_LST", &hdr->target.StartLST, NULL, &status); status=0; 
  ffgky(Fin, TFLOAT, "ION_RM", &hdr->obs.IonRM, NULL, &status); status=0; 
/* Sampling interval (sec) */
  ffgky(Fin, TDOUBLE, "SAMP_INT",&hdr->obs.SampInterval, NULL, &status);status=0; 
/* Output data type String (IQUVLT etc) */
  ffgky(Fin, TSTRING, "OP_STRNG", hdr->obs.OPString, NULL, &status); status=0; 
/* Output scale 1-Arbit; 2-Jy; 3-K;  */
  ffgky(Fin, TINT, "OP_SCALE", &hdr->obs.OPScale, NULL, &status); status=0; 
/* Dynamic DC subtracted? 0/1 */
  ffgky(Fin, TINT, "DYN_DC", &hdr->obs.DynDC, NULL, &status); status=0; 
/* Dynamic DC correction time scale (s) */
  ffgky(Fin, TFLOAT, "DYN_TIME", &hdr->obs.DynTime, NULL, &status); status=0; 
/* No. of bits per sample */
  ffgky(Fin, TINT, "NBITS", &hdr->obs.NBitPerSamp, NULL, &status); status=0; 
/* Interference corrected? 0/1 */
  ffgky(Fin, TINT, "INTERF", &hdr->obs.InterfString, NULL, &status); status=0; 
/* Zap filename */
  ffgky(Fin, TSTRING, "ZAPFILE", hdr->obs.ZapFileName, NULL, &status); status=0; 
/* Number of bins per profile */
  ffgky(Fin, TINT, "NPTSPROF", &hdr->redn.RNBinTimeDump, NULL, &status);status=0;  
/* Number of bins per profile */
  ffgky(Fin, TINT, "NDUMPS", &hdr->redn.RNTimeDumps, NULL, &status); status=0; 
/* Overall Central Sky Frequency */
  ffgky(Fin, TDOUBLE, "FSKYCENT", &hdr->obs.FSkyCent, NULL, &status); status=0; 

  for (i = 0; i < 10; i++)
    {
      ttype[i] = (char *) malloc(20);
      tform[i] = (char *) malloc(20);
      tunit[i] = (char *) malloc(20);
    }
  strcpy(tblname, "BECONFIG");

  ffmnhd(Fin, ASCII_TBL, "BECONFIG", 0, &status);
  ffgnrw(Fin, &nrows, &status);

  indx = 1;
  ffgcv(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.NChan, &anynull, &status);
  indx++;
  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.ChanFreq[i], 
	  &anynull, &status);
    indx++;  
/*     printf("CHANFREQ %d = %f\n",i,hdr->obs.ChanFreq[i]);fflush(stdout); */
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.ChanWidth[i], &anynull, 
	  &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TFLOAT, indx, 1, 1, 1, &fnul, &hdr->obs.ChanTSys[i], &anynull, 
	  &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TFLOAT, indx, 1, 1, 1, &fnul, &hdr->obs.ChanTCal[i], &anynull, 
	  &status);
    indx++;
  }
  ffgcv(Fin, TINT, indx, 1, 1, 1, &fnul, &hdr->asp.SplitFac, &anynull, 
	&status);

  ffmnhd(Fin, ASCII_TBL, "COHDDISP", 0, &status);

  indx = 1;
  ffgcv(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.DM, &anynull, &status); 
  indx++;
/*   printf("DM BEFORE = %f\n",hdr->obs.DM); */
  ffgcv(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.DMMethod, &anynull, &status); 
  indx++;
  ffgcv(Fin, TDOUBLE, indx, 1, 1, 1, &dnul, &hdr->obs.RM, &anynull, &status); 
  indx++;
  ffgcv(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.RMMethod, &anynull, &status); 
  indx++;

  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.ChanChirpLen[i], &anynull,
	  &status);
    indx++;
  }
  for (i=0; i<hdr->obs.NChan; i++) {
    ffgcv(Fin, TINT, indx, 1, 1, 1, &knul, &hdr->obs.ChanOverlap[i], &anynull,
	  &status);
    indx++;
  }

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



  for (i = 0; i < 10; i++)
    {
      free(ttype[i]);
      free(tform[i]);
      free(tunit[i]);
    }

  retval = 0;

  return retval;
}
