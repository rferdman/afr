/* Function to read in parameter table in PSRFITS format file, and output an
   ephemeris par file */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "ASPDefs.h"

int WrtPSRFITSPar(fitsfile *Fin, char *outpar)
{
  int     status=0, anynul;
  int     i_param, n_param, n_f0;
  int     i_col, ncols, colnum, tempo_ver=1, ephver;
  int     tempint;
  long    templng;
  long    nrows;
  float   tempflt;
  double  f0_val;
  double  tempdbl;
  char    junkstr1[128], junkstr2[128];
  char    tempstr1[128], param_ttype[12], param_tform[8];  
  char    ecc_str[8], ecc_val[128];
  char*   temp = &(tempstr1[0]);
  FILE    *Fpar;
  
  /* Open temporary par file for writing */
  if ((Fpar = fopen(outpar, "w")) == NULL){
    printf("Could not open file %s for writing.  Exiting...\n", "temp.par");
    return -1;
  }
  
  /* First see if there exists a PSREPHEM HDU in the FITS file.  If not, 
     then check for PSRPARAM table */

  if(fits_movnam_hdu(Fin, BINARY_TBL, "PSREPHEM", 0, &status)) {
    /* Now try PSRPARAM table */
    status=0;
    if(fits_movnam_hdu(Fin, BINARY_TBL, "PSRPARAM", 0, &status)) {
      fprintf(stderr, "WrtPSRFITSPar ERROR: neither PSREPHEM nor PSRPARAM tables exist!\n");
      return -1;
    }
    else {
      /* Need to read row by row until DM is seen -- and 
	 make sure it's not DM1, etc.*/
      printf("Found PSRPARAM table.\n");

      fits_get_num_rows (Fin, &nrows, &status);
      if(!nrows){
	fprintf(stderr, "WrtPSRFITSPar ERROR: No rows in PSRPARAM table.\n");
	return -1;
      }
      
      /* First get number of rows */
      fits_read_key(Fin, TINT, "NAXIS2", &n_param, NULL, &status); 
      if(nrows!=n_param){
	fprintf(stderr, "WrtPSRFITSPar ERROR: Problem with row numbers in PSRPARAM table.\n");
	return -1;
      }
      /*  will need to read in DM and RM from here as a table column... */
      if(fits_get_colnum(Fin, CASEINSEN, "PARAM", &colnum, &status)){
	fprintf(stderr, "WrtPSRFITSPar ERROR: Problem reading PSRPARAM table.\n");
	return -1;
      }

      //      char* temp = &(tempstr1[0]);
      for (i_param=1; i_param <= n_param; i_param++) {
	//#if 0
	if(fits_read_col(Fin, TSTRING, colnum, i_param, 1, 1, 
			 0, &temp, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Problem reading PSRPARAM table.\n");
	  return -1;
	}
	//#endif
	
	/* Wait on writing out eccentricity until we know tempo version... */
	if(!strncmp(tempstr1, "ECC", 3) || !strncmp(tempstr1, "E ", 2)){
	  sscanf(tempstr1, "%s %s %s", ecc_str, ecc_val, junkstr1);
	}
	/* Need to check EPHVER value and note as a tempo2 par file if needed */
	/* If there is no EPHVER line, it is likely a tempo1 file, and 
	   we have initialized tempo_ver to 1 anyway */
	else if(!strncmp(tempstr1, "EPHVER", 6)){
	  sscanf(tempstr1, "%s %d %s", junkstr1, &ephver, junkstr2);
	  if(ephver < 5){
	    tempo_ver = 1;
	  }
	  else{
	    tempo_ver = 2;	    
	    fprintf(Fpar, "%s\n", tempstr1);
	  }
	}
	else{
	/* Otherwise, write line by line to ascii par file */
	  fprintf(Fpar, "%s\n", tempstr1);
	}

      }
      /* Add in eccentricity line, formatting name according to tempo version*/
      if(tempo_ver == 1){
	fprintf(Fpar, "%-15s %s\n", "E", ecc_val);
      }
      else{
	fprintf(Fpar, "%-15s %s\n", "ECC", ecc_val);	
      }
    }
  }
  else{
    /* If the file contains a PSREPHEM table instead, read and write out 
       par file */
    printf("Found PSREPHEM table.\n");

    /* Set up array of (typical) par file parameters */
    //    char *params[] = {"PSR", "RAJ", "DECJ", "PMRA", "PMDEC", "PX", "F0", 
    //		      "F1", ""}


    /* The PSREPHEM table is a series of columns covering most, if not all, 
       of the parameters that would go into an ephemeris or par file.
       The idea here is to find out the relevant column names with parameter
       values that are nonzero (as best we can) and print the names and 
       values in a temporary par file. */

    /* Use FTGCAL function to find out the ttype and tform of all relevant 
       columns*/
    
    /* First get number of columns is table */
    fits_get_num_cols(Fin, &ncols, &status);
    if(!ncols){
      fprintf(stderr, "WrtPSRFITSPar ERROR: No columns in PSREPHEM table.\n");
      return -1;
    }
    
    /* Read in EPHVER parameter.  If is has a value of >= 5, it is tempo2 format.
       Otherwise, assume it is in tempo1 format. */
    if(fits_get_colnum(Fin, CASEINSEN, "EPHVER", &colnum, &status)){
      fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read PSR_NAME column.\n");
      return -1;
    }

    if(fits_read_col(Fin, TSTRING, colnum, 1, 1, 1, 
		     NULL, &temp, &anynul, &status)) {
      fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read PSR_NAME value.\n");
    }
    printf("tempstr1 = %s\n", tempstr1);
    /* EPHVER can take on weird values if not given in table or is a NULL value
       in the table.  To avoid this, just set its max value to 8, arbitrarily.
       Have set tempo_ver = 1 by default in variable declaration, so if EPHVER
       is NULL or is some weird, unrealistic value, it will default to tempo1 */
    if(tempstr1 != NULL){
      sscanf(tempstr1, "%d", &ephver);
      if(ephver < 8){
	if(ephver < 5)
	  tempo_ver = 1;
	else
	  tempo_ver = 2;
      }    
    }

    /* Know that PSR_NAME corresponds to the PSR parameter of the par file */
    
    if(fits_get_colnum(Fin, CASEINSEN, "PSR_NAME", &colnum, &status)){
      fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read PSR_NAME column.\n");
      return -1;
    }

    if(fits_read_col(Fin, TSTRING, colnum, 1, 1, 1, 
		     NULL, &temp, &anynul, &status)) {
      fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read PSR_NAME value.\n");
    }
    
    /* Write PSR_NAME to file */
    fprintf(Fpar, "%s   %s\n", "PSR", tempstr1);

    /* Begin at next column for the rest */
    f0_val = 0.;
    n_f0 = 0;
    for (i_col=colnum+1; i_col<=ncols; i_col++){
      /* Get column info.  Only care about ttype and tform for now. */
      /* Short form is ftgacl */
      if(fits_get_acolparms(Fin, i_col, &param_ttype[0], NULL, NULL, 
			    &param_tform[0], NULL, NULL, NULL, NULL, &status)){
	fprintf(stderr, "WrtPSRFITSPar ERROR: cannot retrieve information ");
	fprintf(stderr, "for column %d.\n", i_col);
	return -1;
      }
      printf("ttype = %s,  tform = %s\n", param_ttype, param_tform);
      /* If param is eccentricity, make sure it is E for tempo1, and 
	 ECC for tempo2 */
      if(!strcmp(param_ttype, "ECC")){
	if(tempo_ver == 1)
	  strcpy(param_ttype, "E");	
      }

      /* Now read column value, taking into account the data type */
      if(!strncmp(&param_tform[1], "I", 1)){
	if(fits_read_col(Fin, TINT, i_col, 1, 1, 1, NULL, 
			 &tempint, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read (integer) ");
	  fprintf(stderr, "value for column %s.\n", param_ttype);
	}
	/* Now write it out */
	fprintf(Fpar, "%s   %d\n", param_ttype, tempint);
      }
      else if(!strncmp(&param_tform[1], "J", 1)){
	if(fits_read_col(Fin, TLONG, i_col, 1, 1, 1, NULL, 
			 &templng, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read (long integer) ");
	  fprintf(stderr, "value for column %s.\n", param_ttype);
	}
	/* Now write it out unless it is the IF0 column */
	if(strncmp(&param_ttype[0], "IF0", 3)){
	  fprintf(Fpar, "%s   %ld\n", param_ttype, templng);
	}
	else{
	  f0_val += (double)templng;
	  n_f0++; /* increase f0 param counter */
	}
      }
      else if(!strncmp(&param_tform[1], "F", 1)){
	if(fits_read_col(Fin, TFLOAT, i_col, 1, 1, 1, NULL, 
			 &tempflt, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read (float) ");
	  fprintf(stderr, "value for column %s.\n", param_ttype);
	}
	/* Now write it out */
	if(fabsf(tempflt) > FLTTOL) /* Be sure it is non-zero */
	  fprintf(Fpar, "%s   %f\n", param_ttype, tempflt);
      }
      else if(!strncmp(&param_tform[1], "D", 1)){
	if(fits_read_col(Fin, TDOUBLE, i_col, 1, 1, 1, NULL, 
			 &tempdbl, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read (double) ");
	  fprintf(stderr, "value for column %s.\n", param_ttype);
	}
	/* Now write it out unless it is the FF0 column */
	if(strncmp(&param_ttype[0], "FF0", 3)){
	  if(fabs(tempdbl) > FLTTOL)
	    fprintf(Fpar, "%s   %lf\n", param_ttype, tempdbl);
	}
 	else{
	  f0_val += tempdbl;
	  n_f0++; /* increase f0 param counter */
	}
      }
      /* Read anything else as a string */
      else {
	if(fits_read_col(Fin, TSTRING, i_col, 1, 1, 1, NULL, 
			 &temp, &anynul, &status)) {
	  fprintf(stderr, "WrtPSRFITSPar ERROR: Could not read (string) ");
	  fprintf(stderr, "value for column %s.\n", param_ttype);
	}
	fprintf(Fpar, "%s   %s\n", param_ttype, tempstr1);
      }
      
      /* If F0 params are done, write it to file */
      if(!strncmp(&param_ttype[1], "F0", 2) && n_f0 == 2)
	fprintf(Fpar, "%s   %lf\n", "F0", f0_val);	
      
    }
    
  }
  
  fclose(Fpar);
  return 1;
  
}

