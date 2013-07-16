#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPCommon.h"
#include "fitsio.h"

//int ReadASPHdr(struct ASPHdr *, fitsfile *);
//int ReadPSRFITSHdr(struct ASPHdr *, fitsfile *, struct RunVars *);

int ReadHdr(struct ASPHdr *hdr, fitsfile *Fin, struct RunVars *RunMode)
{
  int     hdutype, status=0;


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
/****** DATE GOES HERE ******/ 

/* Backend name ASP/ABPP/WAPP/GUPPI etc */
  fits_read_key(Fin, TSTRING, "BACKEND", hdr->gen.BEName, NULL, &status); 
  status=0; 


  /* If data comes from ASP then do ReadASPData */
  if(!strcmp(hdr->gen.BEName, "xASP")) {
   /* Read in values for header variables */
    if(ReadASPHdr(hdr, Fin) < 0){
      /* Return with error */
      return -1;
    } 
      
  }
  else {  

    /* Otherwise, if data is in PSRFITS format, then do ReadPSRFITSData */
    if(!strcmp(hdr->gen.FitsType, "PSRFITS")) {
      //printf("PSRFITS Data!!\n");

      /* Read in values for header variables */
      if(ReadPSRFITSHdr(hdr, Fin, RunMode) < 0){
	/* Return with error */
	return -1;
      }


    }
    else {
      /* Do not recognize data format! */
      fprintf(stderr, "ReadHdr ERROR: Do not recognize file format.\n");
      return -1;
    }

   }


  return 0;
}
