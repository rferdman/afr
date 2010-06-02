#include <stdio.h>
#include <math.h>
/* #include "fitsio.h" */
#include "ASPCommon.h"

/* int ReadASPData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		fitsfile *, int, int, 
		double **, double **, double **, double **, 
		int **, char **); */

int ReadPSRFITSData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		    fitsfile *, int, int, 
		    double **, double **, double **, double **);

int ReadData(struct ASPHdr *hdr, struct SubHdr *subhdr, 
	     struct RunVars *RunMode, fitsfile *Fin, 
	     int i_dump, int NPtsProf, 
	     double **ASquared, double **BSquared, 
	     double **ReAconjB, double **ImAconjB,
	     int **SampleCount, char **HeadLine)
{


  /* If data comes from ASP then do ReadASPData */
  if(!strcmp(hdr->gen.BEName, "xASP")) {
    //printf("ASP DATA!\n");
    if (ReadASPData(hdr, subhdr, RunMode, Fin, i_dump,
		    hdr->redn.RNBinTimeDump, 
		    ASquared, BSquared, ReAconjB, ImAconjB, 
		    SampleCount, HeadLine) < 0){
      return -1;
    } 
  }
  else {  

    /* Otherwise, if data is in PSRFITS format, then do ReadPSRFITSData */
    if(!strcmp(hdr->gen.FitsType, "PSRFITS")) {
      // printf("PSRFITS Data!!\n");
      if (ReadPSRFITSData(hdr, subhdr, RunMode, Fin, i_dump,
			  hdr->redn.RNBinTimeDump, 
			  ASquared, BSquared, ReAconjB, ImAconjB) < 0){
	/* Return with error */
	return -1;
      }
    }
    else {
      /* Do not recognize data format! */
      fprintf(stderr, "ReadData ERROR: Do not recognize file format.\n");
      return -1;
    }

   }

  return 0;		
}
