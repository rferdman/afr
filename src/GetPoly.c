/* Adapted from GetPoly_master.c in aspire
 * reads info from a polyco.dat file
 * into a Polyco structure.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "polyco.h"
//#include "misc.h"

/* Input:  obs_params needs to be read in before calling this, so
 *         we know the pulsar name we want to find.
 * Output: pc is a array of Polyco structs for every valid polyco set.
 * Returns: number of polyco sets found on success, or -1 on error.
 *          0 means no error, but no valid sets were found.
 */

//int GetPoly(struct asp_params *obs_params, struct Polyco *pc, double ChanFreq) {
int GetPoly(char *polyco_file, char *psr_name, struct Polyco *pc, double ChanFreq, double mjd) {

/*
  Given name,mjd returns coeffs and dm
  Reads polynomial coeffs for current source and time from polyco.dat.
  DM and binary phase info added 23 Apr 90 (new version TZ)
  Puts out two sets of mjdmids for better precision. 16 Nov 1996.
  Recovers Earth-Doppler-shift factor and polyco frequency 13 Feb 1997.
*/

    char name0[15],date0[15],binpha[16];
    char dummy[100];
    char buffer[160];
    /* double aphi0b[30],adphib[30]; */
    double dm0,z40;
    /* double phase; */
    float r;
    long int mjddummy;

    int k,len,jsave;
    int jobs,nblk0,ncoeff0;
    //    double mjd,mjdcheck,mjd1mid,rphase,f0,coeff[16],mjdmid;
    double mjdcheck,mjd1mid,rphase,f0,coeff[16],mjdmid;
    double RefFreq;
    FILE *fp;


    /* Basic check to see that we're on the right day */
    //    mjd = obs_params->imjd + obs_params->fmjd;

    /* Get polyco location */

    /*    char *ptr;
    char polyco_file[256];
    ptr = getenv("HOME");
    if (ptr==NULL) {
        errlog("Coulnd't read HOME env var.\n");
        return(-1);
    } else {
        sprintf(polyco_file, "%s/runtime/polyco.dat", ptr);
	} */

    if ((fp = fopen(polyco_file, "r")) == NULL){
        printf("GetPoly: Could not open %s\n", polyco_file); 
//	errlog("Could not open %s\n", polyco_file); 
       return(-1);
    }
  
    jsave=0;
    
    /* Parse polyco.dat - loop over input */
    while (fgets(buffer,160,fp) != NULL) {
        sscanf(buffer,"%s%s%f%ld%lf%lf%lf",
                name0,date0,&r,&mjddummy,&mjd1mid,&dm0,&z40);
        fscanf(fp,"%lf%lf%d%d%d",&rphase,&f0,&jobs,&nblk0,&ncoeff0);
        fscanf(fp,"%lf%16c",&RefFreq,binpha);
        for(k=0; k<ncoeff0; k++) {
            fscanf(fp,"%s",dummy);
            len = strlen(dummy);
            if (dummy[len-4] == 'D') dummy[len-4] = 'e';
            sscanf(dummy,"%lf",&coeff[k]);
        }
        if (ncoeff0 > MAX_PC_COEFF) {
            printf("GetPoly: ncoeff too big in polyco.dat.\n"); 
//            errlog("ncoeff too big in polyco.dat.\n"); 
            return(-1);
        }
        if (ncoeff0 < MAX_PC_COEFF) {
            for(k=ncoeff0; k<MAX_PC_COEFF; k++) {
                coeff[k] = 0.;
            }
        }
        fgets(buffer,160,fp);  /* skip over line feed */
        mjdmid = mjddummy;
      
        /* just in case its a futmid we want an mjdmid */
        if (mjdmid < 20000) mjdmid += 39126.;
//        if (!strncmp(name0,obs_params->psr_name,10) && (ChanFreq == RefFreq)) {
        if (!strncmp(name0,psr_name,10) && (ChanFreq == RefFreq)) {
            mjdcheck = mjdmid + mjd1mid;
	    //            printf("GETPOLY: mjd = %f, mjdcheck = %f, ChanFreq = %lf, RefFreq = %lf\n",mjd,mjdcheck, ChanFreq, RefFreq); fflush(stdout); 
            if (fabs(mjd-mjdcheck) <= 0.5) {
                pc[jsave].NMinutes = nblk0;
		//                printf("GETPOLY: fabs = %f,   jobs = %d, mjdcheck = %f\n",fabs(mjd-mjdcheck), jobs, mjdcheck); 
		//                printf("GETPOLY: minutes = %d\n",pc[jsave].NMinutes);
                fflush(stdout); 
                pc[jsave].NCoeff = ncoeff0;
                pc[jsave].DM = dm0;
                pc[jsave].EarthZ4 = z40;
                pc[jsave].MjdMidInt = mjdmid;
                pc[jsave].MjdMidFrac = mjd1mid;
                pc[jsave].FRotRef = f0;
                pc[jsave].PhRotRef = rphase;
                pc[jsave].PhRotRef -= floor(pc[jsave].PhRotRef);
                //if ((pc[jsave].PhRotRef < 0.) || 
                //      (pc[jsave].PhRotRef > 1.)) {
                //    /* Phase out of range = do something ??? */
                //}
                for(k=0;k<MAX_PC_COEFF;k++){        
                    pc[jsave].Coeff[k] = coeff[k];
                }
                jsave++;
            }
        } 
    }
    fflush(stdout); 

    /* jsave comes out as the number of polyco sets found for this
     * pulsar.  we'll return this value
     */
    rewind(fp);
    fclose(fp);
    return(jsave);
}


/* Routine to make polyco file, and read it using GetPoly() above, given the pulsar 
   name or par file */
