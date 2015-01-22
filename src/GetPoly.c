/* Adapted from GetPoly_master.c in aspire
 * reads info from a polyco.dat file
 * into a Polyco structure.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
//#include "CmdLine.h"

	 //#define FREQTOL 0.003
#define FREQTOL 0.1

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
	/* Ensure that pulsar names match, and that frequencies match within tolerance */
	//	printf("name0 = %s, psr_name = %s\nfabsf(ChanFreq - RefFreq) = %f, FREQTOL = %f\n", name0, psr_name, fabsf(ChanFreq - RefFreq), FREQTOL);
        if (!strncmp(name0,psr_name,10) && 
	    	(fabsf(ChanFreq - RefFreq) <= FREQTOL) ) {
				mjdcheck = mjdmid + mjd1mid;
	    // printf("GETPOLY: ChanFreq = %lf, RefFreq = %lf\n",ChanFreq, RefFreq); fflush(stdout); 
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


	if(jsave < 1){
		printf("Found a mismatch -- ChanFreq=%.5lf not found...\n",
		ChanFreq);
		printf("Or else an MJD mismatch -- mjd = %lf, mjdcheck = %lf\n", mjd, mjdcheck); 
		printf("Or else a name mismatch -- name0 = %s, psrname = %s\n", name0, psr_name);
		 
	}

    /* jsave comes out as the number of polyco sets found for this
     * pulsar.  we'll return this value
     */
    rewind(fp);
    fclose(fp);
    return(jsave);
}


/* Routine to make polyco file, and read it using GetPoly() above, given the pulsar 
   name or par file -- at the moment this will only be in tempo1 format */

int MakePoly(char* ParFile, struct ASPHdr *Hdr) 
{  
  
  int  i_chan;
  char ParDir[256], tempstr[10];
  //  char ParFile[256];
  char tempo_cmd[256];
  FILE *Fpoly, *Ftempo;
  
  /* Get TEMPO environment variable */
  if(getenv("TEMPO")==NULL){
    fprintf(stderr,"TEMPO environment variable is not set!\n");
    return -1;
  }
  sprintf(ParDir, "%s/tzpar", getenv("TEMPO"));
 
  
  /* Get par file -- Are we using a user-provided par file, or are we looking in the 
     default par directory? */
  
  //  if(Cmd->ParFileP){
    /* Check that input par file exists */
  if(FileExists(ParFile)){
    //      strcpy(ParFile, Cmd->ParFile);
    printf("Creating polycos for PSR %s from input parameter file %s.\n", 
	   Hdr->target.PSRName, ParFile);
  }
  else{
    fprintf(stderr,"Could not open file %s.\n", ParFile);
    return -1;
  }
  //  }
#if 0
  else if(Cmd->PSRNameP){
    printf("Looking for par file in directory %s for pulsar %s, provided by user.\n",
	   ParDir, Cmd->PSRName);
    sprintf(ParFile, "%s/%s.par", ParDir, Cmd->PSRName);
    /* Is there a file in directory that matches Cmd->PSRName.par? */
    if(FileExists(ParFile)){
      printf("Creating polycos from input parameter file %s/%s.par\n",ParDir,Cmd->PSRName);
      /* If yes, take (absolute path of) par file name. */
    }
    else{
      /* If still no, give error and exit. */
      fprintf(stderr, "Could not open file %s.par.  Try explicitly ", Cmd->PSRName);
      fprintf(stderr, "inputting parameter file name by using the -parfile option.\n");
      return -1;
    }
  }
#endif
#if 0
  else{
    printf("Looking for par file in directory %s that matches pulsar name %s ",
	   ParDir, Hdr->target.PSRName);
    printf("from current data file header.\n");
    sprintf(ParFile, "%s/%s.par", ParDir, Hdr->target.PSRName);      
    /* Is there a file in directory that matches Hdr->target.PSRName.par? */
    if(FileExists(ParDir)){
      printf("Creating polycos from input parameter file %s/%s.par\n",
	     ParDir, Hdr->target.PSRName);
      /* If yes, take (absolute path of) par file name. */
    }
    else{
      /* If not, strip first "B" or "J" from name and repeat previous two steps again.*/
      if(!strncmp(Hdr->target.PSRName, "J", 1) || !strncmp(Hdr->target.PSRName, "B", 1)){
	strncpy(tempstr, &Hdr->target.PSRName[1], sizeof(char)*strlen(Hdr->target.PSRName)-1);
	sprintf(ParFile, "%s/%s.par", ParDir, tempstr);
	if(FileExists(ParFile)){
	  printf("Creating polycos from input parameter file %s/%s.par\n",
		 ParDir,tempstr);	  
	  /* If yes, take (absolute path of) par file name. */
	}
	else{
	  fprintf(stderr,"Could not open file %s.par.  Try explicitly ", 
		  tempstr);
	  fprintf(stderr,"inputting parameter file name by using the -parfile option.\n");
	  return -1;	  
	}
      }
      else{
	/* If still no, give error and exit. */
	fprintf(stderr, "Could not open file %s.par.  Try explicitly ", 
		Hdr->target.PSRName);
	fprintf(stderr, "inputting parameter file name by using the -parfile option.\n");
	return -1;
      }
    }
  }
#endif
  
  
  /* If we made it this far, we have a par file we can now use to create a polyco file. */
  
  /* Open final output file and close right away to create it -- sort of 
     like a "touch" command, in effect */
  if((Fpoly=fopen("poly_final.dat", "w"))==NULL){
    fprintf(stderr,"Could not open file poly_final.dat.  \n");
    return -1;
  }
  fclose(Fpoly);
  
  
  /* Loop over number of frequency channels */
  for (i_chan=0; i_chan<Hdr->obs.NChan; i_chan++){
    /* Construct command line based on Header information:  Frequency, 
       start MJD, and get scan length using dump length and number of dumps */
    sprintf(tempo_cmd,
        "tempo -f %s -Zpsr=%s -Zfreq=%lf -Ztobsh=%lf -Zstart=%lf -Zspan=%d -Zsite=%s",
	    /* Channel frequency must be 3 decimal places for tempo1 predictors */
	    ParFile, Hdr->target.PSRName, Hdr->obs.ChanFreq[i_chan], 
	    /* allow an extra 2 hours of polycos */
	    2. + (Hdr->redn.TDump*((double)Hdr->redn.RNTimeDumps)/3600.), 
	    /* Start polyco set to be 1 hour early, (and so end 1 hour later) */
	    Hdr->obs.IMJDStart + ((double)Hdr->obs.StartTime + 
				  floor((double)Hdr->obs.NSubOffs*Hdr->redn.TDump)
				  - 3600.)/86400.,
	    /* default at 15 minutes valid span */
	    1800, Hdr->obs.ObsvtyCode);

    // printf("TDump = %lf,  NDumps = %d\n", Hdr->redn.TDump, Hdr->redn.RNTimeDumps);


    /* Create polyco using tempo command, and append it to final polyco file */
    if(i_chan==1) printf("tempo command:\n%s\n\n\n", tempo_cmd);
    system(tempo_cmd);
#if 0
    if((Ftempo = popen(tempo_cmd, "r"))==NULL) {
      fprintf(stderr, "Could not run tempo correctly\n");
      return -1;
    }
    pclose(Ftempo);
#endif
    /* Now take the output polyco.dat file and append it to the final polyco
       file */
    system("cat polyco.dat >> poly_final.dat");
    if(remove("polyco.dat") < 0) {
      fprintf(stderr,"Could not remove polyco.dat file for frequency %.3lf.",
	      Hdr->obs.ChanFreq[i_chan]);
      return -1;
    }

  }


 return 1;
}




int FileExists(char *testfile)
{
  
  FILE *Ftest;
  
  if((Ftest=fopen(testfile, "r"))==NULL){
    /* File cannot be opened */
    return 0;
  }
  else{
    /* File exists! */
    fclose(Ftest);
    return 1;
  }
  
}

