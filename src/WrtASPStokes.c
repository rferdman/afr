/* ================================================================ */
/* To write the output ASP data in binary tables in FITS format     */
/*                                                                  */
/* INPUTS : Fout    : File Pointer, in "fitsfile" format.           */
/*          nscan   : Current dump                                  */
/*          LSquared, RSquared, ReLconjR, ImLconjR:                 */
/*                    profiles of each quantity, which can be       */
/*                    used to create Stokes parameters              */
/*                                                                  */
/* R. D. Ferdman, 14-January-2004, University of British Columbia   */
/* ================================================================ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPCommon.h"
#include "ASPFitsReader.h"
/* #include "fitsio.h" */

int WrtASPStokes(struct ASPHdr hdr, struct SubHdr subhdr, fitsfile *Fout,
		 int nscan, struct StdProfs *OutputProfs, 
		 struct RunVars *RunMode)
{
  int     status=0, i, j, k, colnum, colindex;
  char    ExtName[4*NCHMAX],*tunit[4*NCHMAX],*tform[4*NCHMAX],*ttype[4*NCHMAX];
  int     retval = -1;
  char    ErrMsg[512];

  double  I[NBINMAX], Q[NBINMAX], U[NBINMAX], V[NBINMAX];

  /* ======================================================================= */
  /* Dump phase and period table                                             */

  /* HERE: add new table for phases calculated due to multiple polycos 
     in each channel*/

  /* table will be freq index, phase, period, and for each dump */
  /* but column number can be freq index, so:  */

  /*
                      
     phase0 period0 
     phase1 period1 
     phase2 period2 
     phase3 period3 
      ...     ...  
      and so on... */

  //printf("NSCAN = %d\n",nscan);fflush(stdout);

  for (i = 0; i < 2; i++){
    ttype[i] = (char *) malloc(64);
    tform[i] = (char *) malloc(64);
    tunit[i] = (char *) malloc(64);
  }
  

  colnum = -1;
  
  sprintf(ttype[++colnum], "REFPHASE");  
  strcpy(tform[colnum],    "F28.12");   /* for ascii table */
  strcpy(tunit[colnum],    "");
  sprintf(ttype[++colnum], "REFPERIOD");  
  strcpy(tform[colnum],    "F20.16");   /* for ascii table */
  strcpy(tunit[colnum],    "s");

  
  strcpy(ExtName,"\0");
  sprintf(ExtName, "DUMPREF%d", nscan); /* One table for each dump */
  if(fits_create_tbl(Fout, ASCII_TBL, 0, ++colnum, ttype, tform, tunit, 
		     ExtName, &status)){
    printf("Could not create first binary table. Error code: %d.  Exiting...\n",
	   status);
    fflush(stdout);
    return -1;
  }

/* Number of seconds since start of obs., taken at middle of dump */
  if(fits_write_key(Fout, TDOUBLE, "MIDSECS", &(subhdr.DumpMiddleSecs), NULL,
		 &status)){
    ffgmsg(ErrMsg);
    printf("Could not write MIDSECS keyword\n");
    printf("  %s\n",ErrMsg);
    fflush(stdout);
    return -1;
  }
  status = 0; 

  colnum = 0;


  if(fits_write_col(Fout, TDOUBLE, ++colnum, 1, 1, (long)hdr.obs.NChan, 
		    &(subhdr.DumpRefPhase[0]), &status)){
    ffgmsg(ErrMsg);
    printf("Could not write REFPHASE column\n");
    printf("  %s\n",ErrMsg);
    fflush(stdout);
    return -1;
  } 
  if(fits_write_col(Fout, TDOUBLE, ++colnum, 1, 1, (long)hdr.obs.NChan,
		    &(subhdr.DumpRefPeriod[0]), &status)){
    ffgmsg(ErrMsg);
    printf("Could not write REFPERIOD column: \n");
    printf("  %s\n",ErrMsg);
    fflush(stdout);
    return -1;
  }
  
  for (i = 0; i < 2; i++){
    free(ttype[i]);
    free(tform[i]);
    free(tunit[i]);
  }

/*======================================================================= */

  /* Place Data table header and table stuff here */
  
  //  printf("WrtASPData says nchan = %d\n",hdr.obs.NChan);fflush(stdout);
  
  for(i=0;i<4*hdr.obs.NChan;i++) {
    ttype[i] = (char*)malloc(64);
    tform[i] = (char*)malloc(64);
    tunit[i] = (char*)malloc(64);
  }

  /* Make Number of tables equal to number of time dumps */

  /* for(i=0;i<SlaveSetup->NDumpsWanted) { */

  /* Set up ttype, tform, and tunit here */
  
  /* Have one set of columns (L^2, R^2, ReLConjR, ImLConjR) 
     for each channel used */
  colnum = -1;
  
  /* Names for each column */
  for(j=0;j<hdr.obs.NChan;j++) {
/* Sort of index each set of column names */
    sprintf(ttype[++colnum], "I%d", j);  
    sprintf(ttype[++colnum], "Q%d", j);
    sprintf(ttype[++colnum], "U%d", j);
    sprintf(ttype[++colnum], "V%d", j);
  }

  
  colnum = -1;

  /* Data format (double precision -- 64 bits) and units (none) 
     for each column */
  for(j=0;j<hdr.obs.NChan;j++) {
   /*  strcpy(tform[++colnum], "F15.10"); */   /* for ascii table */
    strcpy(tform[++colnum], "D");   /* for binary table */
    strcpy(tunit[colnum], "");
    strcpy(tform[++colnum], "D");   /* for binary table */
    strcpy(tunit[colnum], "");
    strcpy(tform[++colnum], "D");   /* for binary table */
    strcpy(tunit[colnum], "");
    strcpy(tform[++colnum], "D");   /* for binary table */
    strcpy(tunit[colnum], "");
  }
  
  strcpy(ExtName,"\0");
  sprintf(ExtName, "STOKES%d", nscan); /* One table for each dump */
  

/*   fits_create_tbl(Fout, ASCII_TBL, 0, 4*SlaveSetup->TotalChansUsed, 
     ttype, tform, tunit, ExtName, &status); */
  if(fits_create_tbl(Fout, BINARY_TBL, 0, 4*hdr.obs.NChan, ttype, tform, tunit,
		     ExtName, &status)){
    printf("Could not create second binary table. Error code: %d.  Exiting...\n",
	   status);
    fflush(stdout);
    return -1;
  }
  /* not sure why nrows is initially set to zero... */


  /* Now write profiles for each channel into table */
  colindex = 1;
  

  //  if(Cmd->OmitfileP) {if(nscan==0){ printf("HELLO AddChans=%d AddDumps=%d!!!\n",RunMode->AddChans,RunMode->AddDumps);fflush(stdout);}}
  for(i=0;i<hdr.obs.NChan;i++) {
    /* If we are not adding *any* dumps OR channels together, check for 
       bad dump/channel combos and convert them to zeroes for later 
       recognition */
    /***** if(Omitfile && !RunMode->AddChans && RunMode->AddDumps==1){
      do_omit=0;
      for(nomit=0;nomit<RunMode->NScanOmit;nomit++){
	if(RunMode->OmitFlag[nscan*hdr.obs.NChan + i]){
	  do_omit=1;
	  break;
	}
      }
      if(do_omit){
	for(k=0;k<hdr.redn.RNBinTimeDump;k++){
	  I[k] = -99999.;
	  Q[k] = -99999.;
	  U[k] = -99999.;
	  V[k] = -99999.;
	}
	do_omit=0;
      }
      else{ 
	for(k=0;k<hdr.redn.RNBinTimeDump;k++){
	  I[k] = (double)(OutputProfs[i].rstds[k]);
	  Q[k] = (double)(OutputProfs[i].rstdq[k]);
	  U[k] = (double)(OutputProfs[i].rstdu[k]);
	  V[k] = (double)(OutputProfs[i].rstdv[k]);
	}
      }
      
    }
    else { ****/
      for(k=0;k<hdr.redn.RNBinTimeDump;k++){
	I[k] = (double)(OutputProfs[i].rstds[k]);
	Q[k] = (double)(OutputProfs[i].rstdq[k]);
	U[k] = (double)(OutputProfs[i].rstdu[k]);
	V[k] = (double)(OutputProfs[i].rstdv[k]);
      }      

      /*****  }  ****/


/* not sure about firstelem, and what it means */
    if(fits_write_col(Fout, TDOUBLE, colindex, 1, 1, 
		      (long)hdr.redn.RNBinTimeDump, &I[0],
		      &status)){
      printf("Could not write column %d into data table.  \n",colindex);
      printf("Error code: %d. Exiting...\n",status);
      fflush(stdout);
      return -2;
    } 
    colindex++; 
    if(fits_write_col(Fout, TDOUBLE, colindex, 1, 1, 
		      (long)hdr.redn.RNBinTimeDump, &Q[0],
		      &status)){
      printf("Could not write column %d into data table.  \n",colindex);
      printf("error code: %d. Exiting...\n",status);
      fflush(stdout);
      return -2;
    } 
    colindex++;
    if(fits_write_col(Fout, TDOUBLE, colindex, 1, 1, 
		      (long)hdr.redn.RNBinTimeDump, &U[0],
		      &status)){
      printf("Could not write column %d into data table.  \n",colindex);
      printf("error code: %d. Exiting...\n",status);
      fflush(stdout);
      return -2;
    }
    colindex++;
    if(fits_write_col(Fout, TDOUBLE, colindex, 1, 1, 
		      (long)hdr.redn.RNBinTimeDump, &V[0],
		      &status)){
      printf("Could not write column %d into data table.  \n",colindex);
      printf("error code: %d. Exiting...\n",status);
      fflush(stdout);
      return -2;
    }
    colindex++;



  }

  for(i=0;i<4*hdr.obs.NChan;i++) {
    free(ttype[i]);
    free(tform[i]);
    free(tunit[i]);
  }
 
  retval = 0;

  return retval;


/*   } */
  
  /* ==================================================================== */




}
