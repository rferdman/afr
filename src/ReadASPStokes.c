#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "ASPCommon.h"

void ReadASPStokes(struct ASPHdr *hdr, struct SubHdr *subhdr,
		   fitsfile *Fin, long NPtsProf, struct StdProfs *Profiles, 
		   int nscan, int Verbose)
{

  int     i, k, status=0, colnum, anynull;
  int     NColumns;

  double  I[NBINMAX], Q[NBINMAX], U[NBINMAX], V[NBINMAX];

  fits_get_num_cols(Fin, &NColumns, &status);status=0;

  if(Verbose)
    printf("NPtsProf = %ld\n\n",NPtsProf);fflush(stdout);
  
  /* get middle dump time  */
  if(!strcmp(hdr->gen.HdrVer,"Ver1.0")) 
    fits_read_key(Fin, TDOUBLE, "DUMPMIDSECS", &(subhdr->DumpMiddleSecs), 
		  NULL, &status); status = 0;


   if(!strcmp(hdr->gen.HdrVer,"Ver1.0.1")) {
     fits_read_key(Fin, TDOUBLE, "MIDSECS",  &(subhdr->DumpMiddleSecs), NULL, 
		   &status); status = 0;
     colnum = 0;
     fits_read_col(Fin, TDOUBLE, ++colnum, 1, 1, (long)hdr->obs.NChan, NULL, 
		   subhdr->DumpRefPhase,  &anynull, &status); 
     fits_read_col(Fin, TDOUBLE, ++colnum, 1, 1, (long)hdr->obs.NChan, NULL, 
		   subhdr->DumpRefPeriod, &anynull, &status); 
     
     if(Verbose){
       printf("Dump %d:  TIME OF DUMP = %lf\n",nscan,subhdr->DumpMiddleSecs); 
       printf("          CHANNEL (MHz)   REF. PHASE   REF. PERIOD (s)\n");
       printf("          -------------   ----------   ---------------\n");
       for(i=0;i<hdr->obs.NChan;i++) 
	 printf("          %13.1lf%13.8lf%20.13lf\n", hdr->obs.ChanFreq[i],
		subhdr->DumpRefPhase[i], subhdr->DumpRefPeriod[i]);       
       fflush(stdout);
     }
     fits_movrel_hdu(Fin, 1, NULL, &status);
   }
   if(Verbose) printf("\n");fflush(stdout);

  colnum = 0;

  for(i=0;i<hdr->obs.NChan;i++) {

    fits_read_col(Fin, TDOUBLE,  ++colnum, 1, 1, NPtsProf, NULL, 
		  &I[0], &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++colnum, 1, 1, NPtsProf, NULL, 
		  &Q[0], &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++colnum, 1, 1, NPtsProf, NULL, 
		  &U[0], &anynull, &status); 
    fits_read_col(Fin, TDOUBLE,  ++colnum, 1, 1, NPtsProf, NULL, 
		  &V[0], &anynull, &status); 

    for (k=0;k<NPtsProf;k++){
      Profiles[i].rstds[k] = (float)(I[k]);
      Profiles[i].rstdq[k] = (float)(Q[k]);
      Profiles[i].rstdu[k] = (float)(U[k]);
      Profiles[i].rstdv[k] = (float)(V[k]);
    }

  }

/* get ref phase and period vectors */
  
  if(!strcmp(hdr->gen.HdrVer,"Ver1.0")){
    printf("FITS header %s\n",hdr->gen.HdrVer);fflush(stdout);
    for(i=0;i<hdr->obs.NChan;i++) {
      fits_read_key(Fin, TDOUBLE, "REFPHASE", &(subhdr->DumpRefPhase[i]),   
		    NULL, &status); status = 0; 
      fits_read_key(Fin, TDOUBLE, "REFPER",   &(subhdr->DumpRefPeriod[i]),  
		    NULL, &status); status = 0;
      /*   *RefPeriod = subhdr->DumpRefPeriod; */
      
      if(Verbose && i==0){
	printf("Dump %d:  TIME OF DUMP = %lf\n",nscan, subhdr->DumpMiddleSecs);
	printf("          REF. PHASE = %lf, REF. PERIOD = %lf s\n", 
	       subhdr->DumpRefPhase[i], subhdr->DumpRefPeriod[i]);fflush(stdout);
      }
    }
  }
  
}
