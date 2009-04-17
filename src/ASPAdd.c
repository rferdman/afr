/* Add up all profiles given in command line -- these are .stokes.fits files
   for now */

/* Only works for single-channel data for now */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"
#include "AddCmdLine.h"

int main(int argc, char **argv)
{

  int             i_file,i_scan,i_chan,i_bin,n_omit,status=0;
  int             NFirstTable, NumHDU, NDumps, TotDumps=0, hdutype;
  // int             OutChans;
  int             spk;
  int             got_bins=0, got_mjd1=0; //, zeroed_outprofs=0;
  long            NPtsProf=0, FirstNPtsProf=0;
  float           Weight=1.0, TotWeight=0.;
  // int             ProfSum=0;
  double          x, ptype;
  double          MJD_first=0., MJD_last=0., MJD_mid;
  double          IMJDMid, MJDSecsMid;
  double          SBase,Srms,Duty,SPeak,FinalMask[NBINMAX];
  double          OutFreq;
  char            ProgName[32];
  char            Outfile[128];
  char            Header[256];  
  struct ASPHdr   *Hdr;
  struct SubHdr   Subhdr;
  struct StdProfs *InProfile, OutProfile;
  struct RunVars  RunMode;
  fitsfile        **Fin;
  FILE            *Fout, *Fcheck;
  Cmdline         *Cmd;


  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  strcpy(ProgName, argv[0]);

  Fin = (fitsfile **)malloc(Cmd->InfileC*sizeof(fitsfile));
  Hdr = (struct ASPHdr *)malloc(Cmd->InfileC*sizeof(struct ASPHdr));

  //    if(!zeroed_outprofs) {
      /*     if(Cmd->SortChansP){
	     OutChans=Hdr[0].obs.NChan;
	     } 
	     else{ */
  // OutChans=1;
  // }
  // OutProfile=(struct StdProfs *)malloc(OutChans*sizeof(struct StdProfs));
  // TotWeight=(float *)malloc(OutChans*sizeof(float));
  
  /* Zero out profiles */
  //     for(i_chan=0;i_chan<OutChans;i_chan++){
  FZero(OutProfile.rstds,NBINMAX);
  FZero(OutProfile.rstdq,NBINMAX);
  FZero(OutProfile.rstdu,NBINMAX);
  FZero(OutProfile.rstdv,NBINMAX);
  // }
  //zeroed_outprofs=1;
      //    }

  /* Create an output file to check omissions if in vebose mode */
  if (Cmd->VerboseP || Cmd->CheckOmitP){
    if((Fcheck = fopen("check_omit.dat","w")) == 0)
       { printf("Cannot open %s. Exiting...\n",Outfile); exit(1); }   
  }


  /* read in all input files and add each to the final profile */

  /* read in all input file headers */
  NPtsProf=0;
  for (i_file=0;i_file<Cmd->InfileC;i_file++){
    n_omit=0;
    status=0;    
    if(fits_open_file(&Fin[i_file], Cmd->Infile[i_file], READONLY, &status)){
      printf("Error opening FITS file %s !!!\n", Cmd->Infile[i_file]);
      exit(1);
    }
    if(ReadASPHdr(&Hdr[i_file], Fin[i_file]) < 0){
      printf("%s> Unable to read Header from file %s.  Exiting...\n",
	     ProgName,Cmd->Infile[i_file]);
      exit(1);
    }

    /* Write file name in verbose-mode omit check file */
    if(Cmd->VerboseP || Cmd->CheckOmitP) 
      fprintf(Fcheck, "\n%s\n",Cmd->Infile[i_file]);


    /* for now just check if all files have same number of channels */
    /* if(Cmd->SortChansP)
        if(i_file>0 && Hdr[i_file].obs.NChan!=Hdr[0].obs.NChan){
        fprintf(stderr,"%s> Different numbers of channels in different files:\n\n",
        ProgName);
        fprintf(stderr,"%s: %d channels,  %s: %d channels\n",
	Cmd->Infile[0],Hdr[0].obs.NChan,
	Cmd->Infile[i_file],Hdr[i_file].obs.NChan);
	} */
    
    /* now find the number of dumps in the file */    
    fits_get_num_hdus(Fin[i_file], &NumHDU, &status);
    if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0")){
      NDumps = NumHDU-3;  /* the "3" is temporary, depending on how 
				   many non-data tables we will be using */
    }
    else if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0.1")){
      NDumps = (NumHDU-3)/2;
    }
    else{
      fprintf(stderr,"%s> Do not recognize FITS file version in header.\n",
	      ProgName);
      fprintf(stderr,"This header is %s. Exiting...\n",Hdr[i_file].gen.HdrVer);
      exit(1);
    }

    printf("File %s:\n",Cmd->Infile[i_file]);
    printf("     Number of channels:  %d\n",Hdr[i_file].obs.NChan) ;
    printf("     Number of dumps:     %d\n",NDumps);

    /* Move to the first data table HDU in the fits file */
    if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0"))
      fits_movnam_hdu(Fin[i_file], BINARY_TBL, "STOKES0", 0, &status);
    else if (!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0.1"))
      fits_movnam_hdu(Fin[i_file], ASCII_TBL, "DUMPREF0", 0, &status);

    /* Get the current HDU number */
    fits_get_hdu_num(Fin[i_file], &NFirstTable);

    /* Set up Profile structure size */


    for(i_scan=0;i_scan<NDumps;i_scan++){
      
      InProfile=(struct StdProfs *)malloc(Hdr[i_file].obs.NChan*
					sizeof(struct StdProfs));
     /* move to next dump's data */
      if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0")){
	fits_movabs_hdu(Fin[i_file], NFirstTable+i_scan, &hdutype, &status); 
      }
      else if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0.1")){
	/* if we've reached the end of the FITS file then increase FileNo */
	fits_movabs_hdu(Fin[i_file],NFirstTable+(i_scan%MAXDUMPS)*2+1,&hdutype,
			&status);
	fits_get_num_rows(Fin[i_file], &NPtsProf, &status);status=0; 
	fits_movrel_hdu(Fin[i_file], -1, NULL, &status);
      }
      
      /* IF not done so, use number of bins from first file to compare to 
	 the rest of the files */
      if(got_bins==0){
	FirstNPtsProf=NPtsProf;  
	got_bins=1;
      }

      /**********  FIX:  SKIP THIS WITHOUT DOING THE OMIT THING ***********/
      /**********         AND DON'T ADD TO NUMBER COUNT  ************/
      if(NPtsProf != FirstNPtsProf) {
	fprintf(stderr,"Warning: Skipping scan %d (%ld bins ",
		i_scan,NPtsProf);
	fprintf(stderr,"vs. %ld bins in others).\n",FirstNPtsProf);
		
	n_omit += Hdr[i_file].obs.NChan;
      }
      else{
	
	/* find NPtsProf */
	ReadASPStokes(&Hdr[i_file], &Subhdr, Fin[i_file], NPtsProf, 
		      InProfile, i_scan, Cmd->VerboseP);
	
	
	/* Add this profile onto running output profile */
	
	for(i_chan=0;i_chan<Hdr[i_file].obs.NChan;i_chan++){
	  
	  /* Bad scans are zeroed so if summ of the profile is zero, it's 
	     not to be used in summation */
	  //	  ProfSum = FSum(&InProfile[i_chan].rstds[0], NPtsProf);
	  //	  ProfSum = 0;
	  // ProfSum = ArrayZero(InProfile[i_chan].rstds, NPtsProf);
	  //	  if(ProfSum != 0.0) { // i.e. good data
	  /* Test that all bins in current profile are not zeroed */
	  if(!ArrayZero(InProfile[i_chan].rstds, NPtsProf)) { // i.e. good data
	    //	  if(InProfile[i_chan].rstds[0] > -99998.) { // i.e. good data
	    
	    /* If first MJD has not been registered, then do so since this 
	       would be the first non-omitted scan */
	    if (got_mjd1==0) {
	      MJD_first = (double)Hdr[i_file].obs.IMJDStart + 
		Subhdr.DumpMiddleSecs/86400.;
	      got_mjd1=1;
	    }
	    if (i_scan==NDumps-1) {
	      /* Just keep overwriting MJD_last every i_file -- that way we 
		 ensure getting the last MJD of the FINAL non-omitted scan 
		 used */
	      MJD_last = (double)Hdr[i_file].obs.IMJDStart + 
		Subhdr.DumpMiddleSecs/86400.;
	    }

	    /* Get SNR for each Profile if we want to use weighting; 
	       otherwise weights will all be 1.0 */
	    if(Cmd->WeightP) {
	      Duty = DutyLookup(Hdr[i_file].target.PSRName);
	      BMask(InProfile[i_chan].rstds,&Hdr[i_file].redn.RNBinTimeDump,
		    &Duty,FinalMask);
	      Baseline(InProfile[i_chan].rstds,FinalMask,
		       &Hdr[i_file].redn.RNBinTimeDump,&SBase,&Srms);
	      SPeak =  FindPeak(InProfile[i_chan].rstds,
				&Hdr[i_file].redn.RNBinTimeDump,&spk);
	      InProfile[i_chan].SNR = SPeak*Srms;
	      //	      Weight = InProfile[i_chan].SNR;
	      Weight = Srms; // which is actually 1/RMS.
	    }
	    // printf("SNR %d = %lf\n",i_chan,InProfile[i_chan].SNR);
	    
	    /* Need to figure out how to organize input channels to match 
	     *  output channels */
	    /* if(Cmd->SortChansP){
	          for(i_bin=0;i_bin<NPtsProf;i_bin++) {
	             OutProfile[i_chan].rstds[i_bin] +=
		     Weight*InProfile[i_chan].rstds[i_bin];
		     OutProfile[i_chan].rstdq[i_bin] +=
		     Weight*InProfile[i_chan].rstdq[i_bin];
		     OutProfile[i_chan].rstdu[i_bin] +=
		     Weight*InProfile[i_chan].rstdu[i_bin];
		     OutProfile[i_chan].rstdv[i_bin] +=
		     Weight*InProfile[i_chan].rstdv[i_bin];
		  }
	       } 
	       else{ */
	    for(i_bin=0;i_bin<NPtsProf;i_bin++) {
	      OutProfile.rstds[i_bin] += 
		Weight*InProfile[i_chan].rstds[i_bin];
	      OutProfile.rstdq[i_bin] += 
		Weight*InProfile[i_chan].rstdq[i_bin];
	      OutProfile.rstdu[i_bin] += 
		Weight*InProfile[i_chan].rstdu[i_bin];
	      OutProfile.rstdv[i_bin] += 
		Weight*InProfile[i_chan].rstdv[i_bin];  
	      //   printf("%f\n",OutProfile[0].rstds[i_bin]);fflush(stdout);
	    }
	    //  }
	    TotWeight += Weight;  // for now keep at zero index
	    /* Print profile weights for each scan, for each channel */
	    if(RunMode.Verbose) {
	      if(i_chan==0) printf("Profile weights -- scan %d: \n   ",i_scan);
	      printf("%6.2f  ",Weight);
	      if(i_chan==Hdr[i_file].obs.NChan-1) printf("\n");fflush(stdout);
	    }
	  }
	  else {
	    n_omit++;
	    if(Cmd->VerboseP || Cmd->CheckOmitP){
	      fprintf(Fcheck, "%6d     %.1lf\n",
		      i_scan,Hdr[i_file].obs.ChanFreq[i_chan]);
	 /* printf("File %d, Dump %d, Channel %d (%lf MHz) were found to be\n",
		   i_file,i_scan,i_chan,Hdr[i_file].obs.ChanFreq[i_chan]);
	    printf("  zeroed and so are not included.\n");fflush(stdout); */
	    }
	  }
	  //    if(i_chan==0) {for(i=0;i<50;i++) printf("%lf  ",OutProfile[0].rstds[i]);printf("\n\n");fflush(stdout);};
	}
	
	/*****************/
	
	/*     }  */
	
	
      } /* else from positive check on NPtsProf */
     free(InProfile);
    }    
    /****** maybe bring this inside the ELSE where entire scans aren't being omitted ******/
    /* if (Cmd->SortChansP) 
      TotDumps += (NDumps - n_omit);
      else */ 
    TotDumps += (NDumps*Hdr[i_file].obs.NChan - n_omit); 
    //free(InProfile);
    
    printf("Reading of file %s complete and successful.\n",
	     Cmd->Infile[i_file]);
    printf("%d scans omitted.\n\n",n_omit);fflush(stdout);
  }

  if(Cmd->VerboseP || Cmd->CheckOmitP) fclose(Fcheck);

  /* Appease the format of the MakePol routine by making up these RunMode 
     structure members */
  strcpy(RunMode.Source,Hdr[0].target.PSRName);
  RunMode.Verbose = Cmd->VerboseP;
  RunMode.FlipPA = 0;
  RunMode.NoBase = Cmd->NoBaseP;
  /* divide out total number of dumps to get the average */
  // for(i_chan=0;i_chan<OutChans;i_chan++){
   
  printf("Totdumps = %d\n",TotDumps);
  fflush(stdout);
  
  for(i_bin=0;i_bin<NPtsProf;i_bin++) {
    OutProfile.rstds[i_bin] /= TotWeight;
    OutProfile.rstdq[i_bin] /= TotWeight;
    OutProfile.rstdu[i_bin] /= TotWeight;
    OutProfile.rstdv[i_bin] /= TotWeight;
  }
  
  MakePol(&RunMode, (int)NPtsProf, &OutProfile);
  /* Open file for writing */
  sprintf(Outfile,"AddProf.out");
  
  /* now write the output ascii added profile */
  if ((Fout = fopen(Outfile,"w")) == 0)
    { printf("Cannot open %s. Exiting...\n",Outfile); exit(1); }
  
  /* take average MJD of first to last scan */
  MJD_mid = (MJD_first + MJD_last)/2.;
  IMJDMid = floor(MJD_mid);
  MJDSecsMid = (MJD_mid - IMJDMid)*86400.;
  
  printf("MJD_mid = %lf, IMJDMid = %lf, MJDSecsMid = %lf\n",MJD_mid,IMJDMid,MJDSecsMid);fflush(stdout); 

  /* to choose a channel to put in the header for now, ust use the average 
     of the first datafile's channels */ 
  OutFreq=0.;
  for(i_chan=0; i_chan<Hdr[0].obs.NChan; i_chan++){
    OutFreq += Hdr[0].obs.ChanFreq[i_chan];
  }
  OutFreq /= Hdr[0].obs.NChan;

  /* Create and print header line for output file(s) */
  sprintf(Header,"# %.1f %.7f %.10f %ld %.3f %.3f %d %s %d %9s %.10f",
	  IMJDMid, 
	  MJDSecsMid, 
	  //	  Subhdr.DumpRefPeriod[i_chan],
	  0.,
	  (long)1, OutFreq,
	  //	    Hdr[0].obs.DM, Hdr[i_file].redn.RNBinTimeDump,
	  Hdr[0].obs.DM, (int)NPtsProf,
	  Hdr[0].obs.ObsvtyCode, 1, Hdr[0].target.PSRName, 
	  0.);             
  //	  Subhdr.DumpRefPhase[i_chan]);             
  fprintf(Fout,"%s\n",Header);
  
  for(i_bin=0;i_bin<NPtsProf;i_bin++) {
    /* see how strong the linear polarization is */
    x = OutProfile.stdlin[i_bin]*OutProfile.Srms; 
    ptype = 43.1;
    if (x > 1.) ptype=43.2;
    if (x > 2.) ptype=43.3;
    if (x > 3.) ptype=43.4;
    if (x > 4.) ptype=43.5;
    if (x > 5.) ptype=43.6;
    fprintf(Fout,"%5d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%6.1f\n",i_bin,
	    OutProfile.rstds[i_bin],OutProfile.rstdq[i_bin],
	    OutProfile.rstdu[i_bin],
	    OutProfile.rstdv[i_bin],
	    /* phi in degrees */
	    OutProfile.stdlin[i_bin],
	    OutProfile.stdphi[i_bin]*180.0/TWOPI, 
	    OutProfile.stdphierr[i_bin]*180.0/TWOPI,ptype);
  }
  printf("Created output file %s\n",Outfile);
  fclose(Fout);
  // }
  
  /* Write all this to file */
  
  
  printf("\nCompleted successfully.\n\n");fflush(stdout);
  
  exit(0);

}

