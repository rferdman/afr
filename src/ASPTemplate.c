/* Make a standard profile with the option of giving an ascii profile to 
   which to add *.stokes.fits files 

   -- RDF, 3 April, 2009 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"
#include "TemplateCmdLine.h"

int main(int argc, char **argv)
{

  int             i_file,i_dump,i_chan,i_bin,n_omit,status=0;
  int             NFirstTable, NumHDU, NDumps, TotDumps=0, hdutype;
  int             spk;
  int             got_bins=0, got_mjd1=0, good_freqs=1, bad_freqs;
  int             bin[NBINMAX], NStdBins;
  int             ngood;

  long            NPtsProf=0, FirstNPtsProf=0;

  float           Weight=1.0, TotWeight=0.;
  float           profs[NBINMAX],amps[NBINMAX], phas[NBINMAX];
  float           Shift,EShift,SNR,ESNR,b,errb;
  float           CutThresh=0.; // degree of Fourier-space noise cutting

  double          x, ptype;
  double          MJD_first=0., MJD_last=0., MJD_mid;
  double          IMJDMid, MJDSecsMid;
  double          SBase,Srms,Duty,SPeak,FinalMask[NBINMAX];
  double          FreqRange[2], OutFreq;
  double          ByAngle;

  char            ProgName[32];
  char            Outfile[128];
  char            Header[256];  
  char            Headerline[256];
  char            PSRName[16];

  struct ASPHdr   *Hdr;
  struct SubHdr   Subhdr;
  struct StdProfs *InProfile, OutProfile;
  struct StdProfs StdProfile;
  struct StdProfs TempProf; 
  struct RunVars  RunMode, TempRunMode;

  fitsfile        **Fin;
  FILE            *Fout, *Fcheck;
  Cmdline         *Cmd;

  
  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  strcpy(ProgName, argv[0]);


  /* Make sure that at least one of standard profile or new input profiles
     have been provided */
  if (!(Cmd->TemplateP || Cmd->InfileP)) {
    printf("Must select one of standard template profile or input ");
    printf("*.fits profiles, or both.  Exiting...\n");
    exit(1);
  }

  /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  RunMode.Verbose = Cmd->VerboseP;
  RunMode.FlipPA = 0;
  RunMode.NoBase = Cmd->NoBaseP;
  //    if(!zeroed_outprofs) {
  /*     if(Cmd->SortChansP){
	 OutChans=Hdr[0].obs.NChan;
	 } 
	 else{ */

  /* Check whteher user has requested a frequency range.  If so, set FreqRange array.
     If not, make FreqRange default to large boundaries. */
  FreqRange[0] = 0.;
  FreqRange[1] = 9999999.;
  if(Cmd->FreqRangeP){
    if(Cmd->FreqRange[0] < Cmd->FreqRange[1]) {
      FreqRange[0] = Cmd->FreqRange[0];
      FreqRange[1] = Cmd->FreqRange[1];
    }
    else {
      FreqRange[0] = Cmd->FreqRange[1];
      FreqRange[1] = Cmd->FreqRange[0];
    }
    printf("Will only add together data within the frequency range %.2lf to %.2lf MHz.\n",
	   FreqRange[0], FreqRange[1]);
  }


  if(Cmd->NoiseCutP) {
    if(Cmd->NoiseCutC == 0) CutThresh = 0.; // minimum of minimizing function
    else if (Cmd->NoiseCut < 0) CutThresh = -1.;
    else if (Cmd->NoiseCut > 1.0) {
      printf("Values of -noisecut argument must be less than 1.0.  ");
      printf("Exiting...\n");
      exit(3);
    }
    else CutThresh = Cmd->NoiseCut;
  }


    

  // }
  // OutProfile=(struct StdProfs *)malloc(OutChans*sizeof(struct StdProfs));
  // TotWeight=(float *)malloc(OutChans*sizeof(float));
  
  /* Zero out profiles */
  //     for(i_chan=0;i_chan<OutChans;i_chan++){
  FZero(OutProfile.rstds,NBINMAX);
  FZero(OutProfile.rstdq,NBINMAX);
  FZero(OutProfile.rstdu,NBINMAX);
  FZero(OutProfile.rstdv,NBINMAX);

  
  /* Create an output file to check omissions if in vebose mode */
  if (Cmd->VerboseP || Cmd->CheckOmitP){
    if((Fcheck = fopen("check_omit.dat","a")) == 0)
      { printf("Cannot open check_omit.dat. Exiting...\n"); exit(1); }   
  }


  /* If user uses the -template option then we will be scaling all profiles 
     using fftfit before adding them to the total profile  */
  // if (Cmd->TemplateP){
  if(Cmd->RotateP) 
    printf("Will rotate input profiles to match template profile.\n\n");
  

  /* Now, if the user has provided an ascii profile to scale things by, 
     read it in, make sure baseline is subtracted off, 
     and get Fourier compnents */
  if (Cmd->TemplateP){
    printf("Will scale input and output profiles, if requested, by the ");
    printf("found in %s.\n\n", Cmd->Template);
    if (!Cmd->PulsarP){
      printf("Must identify pulsar name with the -psr command when using ");
      printf("input standard profile.  Exiting...\n");
      exit(2);
    }
    
    if (ReadASPAsc(Cmd->Template, Headerline, bin, &StdProfile, &NStdBins) < 0) {
      printf("Error in reading file %s.\n",Cmd->Template);
      fflush(stdout);
      exit(1);
    }
    strcpy(PSRName,Cmd->Pulsar);
    strcpy(RunMode.Source,Cmd->Pulsar);
    RemoveBase(&RunMode, NStdBins, &StdProfile);
  }
    /* If Template was provided, then set up fourier space components 
       for fftfit'ing later */
  if (Cmd->TemplateP){
    cprofc(StdProfile.rstds,NStdBins,
	   StdProfile.stdamps,StdProfile.stdphas);
    
    memcpy(amps,StdProfile.stdamps,sizeof(float)*NBINMAX);
    memcpy(phas,StdProfile.stdphas,sizeof(float)*NBINMAX);     
    
  }

  if (!Cmd->InfileP){
    /* First, make sure that we are doing a noise cut, otherwise there is 
       no point to copying it over! */
    if (!Cmd->NoiseCutP) {
      printf("If only using standard profile, then must perform a noise cut ");
      printf("with -noisecut option.  Exiting....\n");
      exit(2);
    }
    /* If no Infiles given, then just copy StdProfile to OutProfile */
    memcpy(&OutProfile, &StdProfile, sizeof(struct StdProfs));
    memcpy(Header, Headerline, 256); // both are 256-character strings
    NPtsProf = NStdBins;
    /* Appease the format of the MakePol routine by making up these RunMode 
       structure members */
    MakePol(&RunMode, (int)NPtsProf, &OutProfile);
  }

#if 0      
  else {

    /* If Template was provided, then set up fourier space components 
       for fftfit'ing later */
    if (Cmd->TemplateP){
      cprofc(StdProfile.rstds,NStdBins,
	     StdProfile.stdamps,StdProfile.stdphas);
      
      memcpy(amps,StdProfile.stdamps,sizeof(float)*NBINMAX);
      memcpy(phas,StdProfile.stdphas,sizeof(float)*NBINMAX);     

      /* Now initialize the Output profile, starting with this input
	 standard profile, first taking into account the weight */
      if(!Cmd->NoWeightP) {
	Duty = DutyLookup(PSRName);
	BMask(StdProfile.rstds, &NStdBins, &Duty,FinalMask);
	Baseline(StdProfile.rstds,FinalMask, &NStdBins,&SBase,&Srms);
	SPeak =  FindPeak(StdProfile.rstds, &NStdBins, &spk); 
	StdProfile.SNR = (SPeak-SBase)*Srms;
	Weight = StdProfile.SNR;
	//Weight = Srms; // which is actually 1/RMS.
	//printf("Standard profile weight = %.3lf\n", Weight);
      }
      else {
	Weight = 1.0;
      }
      /* If we want to pre-scale the current input profile  */
      if (Cmd->ScaleBeforeP) {
	if(Cmd->ScaleBeforeC == 0){		    
	  Zero2One(&StdProfile, NStdBins, PSRName);	  
	}
	/* Scale profile to input template before adding to output profile */
	else{
	  
	}
      }
      	      
      for(i_bin=0;i_bin<NStdBins;i_bin++) {
	OutProfile.rstds[i_bin] += 
	  Weight*StdProfile.rstds[i_bin];
	OutProfile.rstdq[i_bin] += 
	  Weight*StdProfile.rstdq[i_bin];
	OutProfile.rstdu[i_bin] += 
	  Weight*StdProfile.rstdu[i_bin];
	OutProfile.rstdv[i_bin] += 
	  Weight*StdProfile.rstdv[i_bin];  
	//   printf("%f\n",OutProfile[0].rstds[i_bin]);fflush(stdout);
      }
      TotWeight += Weight;
      
    }
    // }
#endif



    /* If infiles not given then just copy input standard profile to output
       profile */

    /* read in all input files and add each to the final profile */
    Fin = (fitsfile **)malloc(Cmd->InfileC*sizeof(fitsfile));
    Hdr = (struct ASPHdr *)malloc(Cmd->InfileC*sizeof(struct ASPHdr));
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


      /* Do check of Frequency range here.  If NONE of profile is within range, then there
	 is no point in continuing to read this file. For now ALL frequencies MUST be within
	 range, so user must choose range wisely.. */
      if(Cmd->VerboseP){
	printf("\nTEST:  ChanFreq[0]  = %.3lf,  ChanFreq[last] = %.3lf\n",
	       Hdr[i_file].obs.ChanFreq[0], 
	       Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);
	printf("TEST:  FreqRange[0] = %.3lf,  FreqRange[1]   = %.3lf\n\n",
	       FreqRange[0], FreqRange[1]);
	fflush(stdout);
      }
      
      if(Hdr[i_file].obs.ChanFreq[0] <= Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]){
	printf("\n\nFrequency range = %.3lf --> %.3lf\n", 
	       Hdr[i_file].obs.ChanFreq[0],
	       Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);
	bad_freqs = ( (FreqRange[0] < Hdr[i_file].obs.ChanFreq[0]) && 
		      (FreqRange[1] < Hdr[i_file].obs.ChanFreq[0]) ) ||
	  ( (FreqRange[0] > Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]) && 
	    (FreqRange[1] > Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]) );
	
	//	good_freqs = (FreqRange[0] >= Hdr[i_file].obs.ChanFreq[0] ||
	//	      FreqRange[1] <= Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);
      }
      else{
	printf("\n\nFrequency range = %.3lf --> %.3lf\n",
	       Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1],
	       Hdr[i_file].obs.ChanFreq[0]);
	bad_freqs = ( (FreqRange[0] > Hdr[i_file].obs.ChanFreq[0]) && 
		      (FreqRange[1] > Hdr[i_file].obs.ChanFreq[0]) ) ||
	  ( (FreqRange[0] < Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]) && 
	    (FreqRange[1] < Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]) );
	//good_freqs = (FreqRange[1] <= Hdr[i_file].obs.ChanFreq[0] ||
	//	      FreqRange[0] >= Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);

      }
      //      printf("\n\nFrequency range = %.3lf --> %.3lf\n", Hdr[i_file].obs.ChanFreq[0],
      //	     Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);
      good_freqs = !bad_freqs;
      printf("Frequency restrict range: %.3lf --> %.3lf\n", FreqRange[0], FreqRange[1]);
      printf("good_freqs = %d\n\n", good_freqs);

      /* Set Source Name here */
      if(i_file==0) {
	strcpy(RunMode.Source,Hdr[0].target.PSRName);
	strcpy(PSRName,Hdr[0].target.PSRName);
      }	


      /* Write file name in verbose-mode omit check file */
      /* if(Cmd->VerboseP || Cmd->CheckOmitP) 
	 fprintf(Fcheck, "\n%s\n",Cmd->Infile[i_file]); */


      /* for now just check if all files have same number of channels */
      /* if(Cmd->SortChansP)
	 if(i_file>0 && Hdr[i_file].obs.NChan!=Hdr[0].obs.NChan){
	 fprintf(stderr,"%s> Different numbers of channels in different files:\n\n",
	 ProgName);
	 fprintf(stderr,"%s: %d channels,  %s: %d channels\n",
	 Cmd->Infile[0],Hdr[0].obs.NChan,
	 Cmd->Infile[i_file],Hdr[i_file].obs.NChan);
	 } */
    
      /* If observing frequencies are within restricted range, then go ahead */
      if(good_freqs) {

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
	
	
	for(i_dump=0;i_dump<NDumps;i_dump++){
	  
	  /* Set up Profile structure size */
	  InProfile=(struct StdProfs *)malloc(Hdr[i_file].obs.NChan*
					      sizeof(struct StdProfs));
	  
	  /* move to next dump's data */
	  if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0")){
	    fits_movabs_hdu(Fin[i_file], NFirstTable+i_dump, &hdutype, &status); 
	  }
	  else if(!strcmp(Hdr[i_file].gen.HdrVer,"Ver1.0.1")){
	    /* if we've reached the end of the FITS file then increase FileNo */
	    fits_movabs_hdu(Fin[i_file],NFirstTable+(i_dump%MAXDUMPS)*2+1,
			    &hdutype,  &status);
	    fits_get_num_rows(Fin[i_file], &NPtsProf, &status);status=0; 
	    fits_movrel_hdu(Fin[i_file], -1, NULL, &status);
	  }
	  
	  /* IF not done so, use number of bins from first file to compare to 
	     the rest of the files */
	  if(got_bins==0){
	    FirstNPtsProf=NPtsProf;  
	    /* Set number of bins in first file's profiles to be StdBins if we 
	       don't input a standard profile */
	    if(!Cmd->TemplateP) NStdBins=NPtsProf;
	    got_bins=1;
	  }
	  
	  /**********  FIX:  SKIP THIS WITHOUT DOING THE OMIT THING ***********/
	  /**********         AND DON'T ADD TO NUMBER COUNT  ************/
	  /* if(NPtsProf != NStdBins) {
	    fprintf(stderr,"Warning: Skipping scan %d (%ld bins ",
		    i_dump,NPtsProf);
	    fprintf(stderr,"vs. %d bins in others).\n",NStdBins);
	    
	    n_omit += Hdr[i_file].obs.NChan;
	    } 
	  else{ */
	    
	  /* find NPtsProf */
	  ReadASPStokes(&Hdr[i_file], &Subhdr, Fin[i_file], NPtsProf, 
			InProfile, i_dump, Cmd->VerboseP);
	  
	  /* Add this profile onto running output profile */
	  
	  for(i_chan=0;i_chan<Hdr[i_file].obs.NChan;i_chan++){
	    
	    /* First filter -- make sure that the current channel falls within desired 
	       frequency range restriction */
	    if(Hdr[i_file].obs.ChanFreq[i_chan] >= FreqRange[0] &&
	       Hdr[i_file].obs.ChanFreq[i_chan] <= FreqRange[1]) { 
	      
	      
	      /* Bad scans are zeroed so if summ of the profile is zero, it's 
		 not to be used in summation */
	      if(!ArrayZero(InProfile[i_chan].rstds, NPtsProf)) { // i.e. good data
		if(!Cmd->NoBaseP) RemoveBase(&RunMode, NPtsProf, &InProfile[i_chan]);
		
		/* If first MJD has not been registered, then do so since this 
		   would be the first non-omitted scan */
		if (got_mjd1==0) {
		  MJD_first = (double)Hdr[i_file].obs.IMJDStart + 
		    Subhdr.DumpMiddleSecs/86400.;
		  got_mjd1=1;
		}
		if (i_dump==NDumps-1) {
		  /* Just keep overwriting MJD_last every i_file -- that way we 
		     ensure getting the last MJD of the FINAL non-omitted scan 
		     used */
		  MJD_last = (double)Hdr[i_file].obs.IMJDStart + 
		    Subhdr.DumpMiddleSecs/86400.;
		}
		
		// printf("SNR %d = %lf\n",i_chan,InProfile[i_chan].SNR);
		
		/* Now, if std profile is provided, then scale current scan to it 
		   using fftfit.  If std profile is NOT provided, then let first 
		   file accumulate and copy that into StdProfile, scaling 
		   subsequent profiles to that from then on. */
		
		if(Cmd->TemplateP || (!Cmd->TemplateP && i_file > 0)) {
		  //	      if (Cmd->TemplateP || i_file > 0){
		  /* Scale by scale factor unless it is the first file and no 
		     standard profile is given */
		  /* Run fftfit */
		  if (Cmd->ScaleP || Cmd->RotateP){
		    memcpy(profs,InProfile[i_chan].rstds,sizeof(float)*NBINMAX);
		    /* Now fftfit to find shift required in second profile */
		    fftfit_(profs,&amps[1],&phas[1],
			    &NStdBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);
		  }
		  
		  if (Cmd->ScaleP) {
		    //printf("Scaling according to standard profile...\n");
		    printf("scale factor = %f\n\n", b);fflush(stdout);
		    for(i_bin=0; i_bin<NPtsProf; i_bin++){
		      InProfile[i_chan].rstds[i_bin] /= b;
		      InProfile[i_chan].rstdq[i_bin] /= b;
		      InProfile[i_chan].rstdu[i_bin] /= b;
		      InProfile[i_chan].rstdv[i_bin] /= b;
		    } 
		  }
		  
		  if(Cmd->RotateP) {
		    RunMode.Verbose = Cmd->VerboseP;
		    RunMode.NBins = NPtsProf;
		    
		    ByAngle = (double)(-Shift/NPtsProf*TWOPI);
		    if(Cmd->VerboseP){
		      printf("Rotating dump %d, channel %d (%lf.1 MHz)\n", 
			     i_dump, i_chan, Hdr[i_file].obs.ChanFreq[i_chan]);
		      printf("Scale factor: %f +- %f \n",b,errb);
		    }
		    RotateProf(&RunMode, &InProfile[i_chan], ByAngle);
		  }
		}
		/* Get SNR for each Profile if we want to use weighting; 
		   otherwise weights will all be 1.0 */
		if(!Cmd->NoWeightP) {
		  Duty = DutyLookup(PSRName);
		  BMask(InProfile[i_chan].rstds,&Hdr[i_file].redn.RNBinTimeDump,
			&Duty,FinalMask);
		  Baseline(InProfile[i_chan].rstds,FinalMask,
			   &Hdr[i_file].redn.RNBinTimeDump,&SBase,&Srms);
		  SPeak =  FindPeak(InProfile[i_chan].rstds,
				    &Hdr[i_file].redn.RNBinTimeDump,&spk);
		  InProfile[i_chan].SNR = (SPeak-SBase)*Srms;
		  Weight = InProfile[i_chan].SNR;
		  //Weight = Srms; // which is actually 1/RMS.
		  //printf("Weight = %.2f\n", Weight);
		}
		/* If we want to pre-scale the current input profile  */
		if (Cmd->ScaleBeforeP) {
		  /* Scale profile to input template before adding to output profile */
		  if(Cmd->TemplateP){	
		    /* Check if input profile has same bins as template.  If not,
		       bin downinput profile (or template) for the purpose of
		       calculating scale factor, then apply it to original input 
		       profile */
		    if(NPtsProf > NStdBins){
		      printf("Template profile has fewer bins that input profile(s).  ");
		      printf("Binning down data to match template before scaling.\n\n");
		      TempRunMode.BinDown = 1;
		      TempRunMode.NBins = NPtsProf;
		      TempRunMode.NBinsOut = NStdBins;
		      TempRunMode.Verbose = Cmd->VerboseP;
		      BinDown(&TempRunMode, &InProfile[i_chan], &TempProf);
		      memcpy(profs,TempProf.rstds,sizeof(float)*NBINMAX);		      
		    }
		    else if (NPtsProf < NStdBins) {
		      fprintf(stderr, "Input template has more bins that input ");
		      fprintf(stderr, "profile(s).  Please use binned down version");
		      fprintf(stderr, "of template.  Exiting.\n");
		      exit(3);
		    }
		    /* else, they are equal, so no binning down needed! */
		    else{
		      memcpy(profs, InProfile[i_chan].rstds, sizeof(float)*NBINMAX);		      
		    }
		    if (i_dump==0 && i_chan==0) 
		      printf("Scaling input profiles to match template...\n\n");
		    /* Now fftfit to find shift required in second profile */
		    fftfit_(profs,&amps[1],&phas[1],
			    &NStdBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);
		    for(i_bin=0; i_bin<NPtsProf; i_bin++){
		      /* Divide by scale factor, and by the ratio of bins, so that power 
		       is properly distributed back in non-binned-down profile */
		      InProfile[i_chan].rstds[i_bin] /= (b/(double)(NPtsProf/NStdBins));
		      InProfile[i_chan].rstdq[i_bin] /= (b/(double)(NPtsProf/NStdBins));
		      InProfile[i_chan].rstdu[i_bin] /= (b/(double)(NPtsProf/NStdBins));
		      InProfile[i_chan].rstdv[i_bin] /= (b/(double)(NPtsProf/NStdBins));
		    } 
		  }
		  else{
		    /* If we don't want to scale against a template, simply normalize
		       profile to run between [0,1] */
		    Zero2One(&InProfile[i_chan], NPtsProf, PSRName);		      
		  }
		}
		
		/* Now add onto accumulating output profile */
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
		if(Cmd->VerboseP) {
		  // if(i_chan==0) 
		  printf("Profile weights -- scan %d, chan %d: \n   ",i_dump, i_chan);
		  printf("%6.2f \n ",Weight);
		  if(i_chan==Hdr[i_file].obs.NChan-1) printf("\n");fflush(stdout);
		}
	      }
	      else {
		n_omit++;
		if(Cmd->VerboseP || Cmd->CheckOmitP){
		  fprintf(Fcheck, "%s -- %6d     %.1lf\n",
			  Cmd->Infile[i_file], i_dump,Hdr[i_file].obs.ChanFreq[i_chan]);
		}
	      }
	      
	    } /* if (within frequency range) */
	    else {
	      if(Cmd->VerboseP){
		printf("Scan %d, channel %d (%.2lf MHz) omitted, since it does not fall within ",
		       i_dump, i_chan, Hdr[i_file].obs.ChanFreq[i_chan]);
		printf("specified frequency range.\n");
	      }
	    }
	    
	  }
	  
	  /*****************/
	  
	  /*     }  */
	  
	  
	  //} /* else from positive check on NPtsProf */
	  free(InProfile);
	  
	  /* Now if we *aren't* using a user-supplied standard profile, copy 
	     the current output profile's worth of added scans to the 
	     standard profile for scaling purposes from now on */
	  
	  /* To ensure that this starts with the first file for which all 
	     scans are not omitted, and only for when standard profile is not 
	     given, and do once every dump, not once every channel */
	  if (got_mjd1 && !Cmd->TemplateP) {
	    
	    //printf("Copying over std prof...\n\n");fflush(stdout);
	    for(i_bin=0;i_bin<NPtsProf;i_bin++) {
	      StdProfile.rstds[i_bin] = OutProfile.rstds[i_bin];///TotWeight;
	      StdProfile.rstdq[i_bin] = OutProfile.rstdq[i_bin];///TotWeight;
	      StdProfile.rstdu[i_bin] = OutProfile.rstdu[i_bin];///TotWeight;
	      StdProfile.rstdv[i_bin] = OutProfile.rstdv[i_bin];///TotWeight;
	    }
	    //memcpy(&StdProfile, &OutProfile, sizeof(struct StdProfs));
	    
	    /* Now prepare it for fftfitting next round */
	    cprofc(StdProfile.rstds,NStdBins,
		   StdProfile.stdamps,StdProfile.stdphas);
	    
	    memcpy(amps,StdProfile.stdamps,sizeof(float)*NBINMAX);
	    memcpy(phas,StdProfile.stdphas,sizeof(float)*NBINMAX);     
	  }
	}
	
	/****** maybe bring this inside the ELSE where entire scans aren't being omitted ******/
	/* if (Cmd->SortChansP) 
	   TotDumps += (NDumps - n_omit);
	   else */ 
	TotDumps += (NDumps*Hdr[i_file].obs.NChan - n_omit); 
	
	printf("Reading of file %s complete and successful.\n",
	       Cmd->Infile[i_file]);
	printf("%d scans omitted.\n\n",n_omit);fflush(stdout);
	
	
      } /* if(good_freq) */
      else{
	printf("File %s contains data outside frequency range restrictions ",
	       Cmd->Infile[i_file]);
	printf("(%.2lf to %.2lf).  Skipping this file in constructing template profile.\n",
	       Hdr[i_file].obs.ChanFreq[0], Hdr[i_file].obs.ChanFreq[Hdr[i_file].obs.NChan-1]);
      }
      
    }
    
    if(Cmd->VerboseP || Cmd->CheckOmitP) fclose(Fcheck);
    
    /* Appease the format of the MakePol routine by making up these RunMode 
       structure members */
    /* divide out total number of dumps to get the average */
    // for(i_chan=0;i_chan<OutChans;i_chan++){
    
    printf("Totdumps = %d\nTotWeight = %f\n", TotDumps, TotWeight);
    fflush(stdout);
    
    if (Cmd->NormalizeP){
      for(i_bin=0;i_bin<NStdBins;i_bin++) {
	OutProfile.rstds[i_bin] /= TotWeight;
	OutProfile.rstdq[i_bin] /= TotWeight;
	OutProfile.rstdu[i_bin] /= TotWeight;
	OutProfile.rstdv[i_bin] /= TotWeight;
      } 
    }
  
    MakePol(&RunMode, (int)NPtsProf, &OutProfile);
    /* If we want to scale the output profile  */
    if (Cmd->ScaleAfterP) {
      /* If we don't want to scale against a template, simply normalize
	 profile to run between [0,1] */
      if(Cmd->TemplateP){		    
	if(NPtsProf > NStdBins){
	  printf("Template profile has fewer bins that final output profile.  ");
	  printf("Binning down data to match template before scaling.\n\n");
	  TempRunMode.BinDown = 1;
	  TempRunMode.NBins = NPtsProf;
	  TempRunMode.NBinsOut = NStdBins;
	  TempRunMode.Verbose = Cmd->VerboseP;
	  BinDown(&TempRunMode, &OutProfile, &TempProf);
	  memcpy(profs, TempProf.rstds, sizeof(float)*NBINMAX);		      
	}
	else if (NPtsProf < NStdBins) {
	  fprintf(stderr, "Input template has more bins that output ");
	  fprintf(stderr, "profile.  Please use binned down version");
	  fprintf(stderr, "of template.  Exiting.\n");
	  exit(3);
	}
	/* else, they are equal, so no binning down needed! */
	else{
	  memcpy(profs,OutProfile.rstds,sizeof(float)*NBINMAX);
	}
	/* Now fftfit to find shift required in second profile */
	printf("Scaling output profile to match template...\n\n");
	fftfit_(profs,&amps[1],&phas[1],
		&NStdBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);
	for(i_bin=0; i_bin<NPtsProf; i_bin++){
	  OutProfile.rstds[i_bin] /= (b/(double)(NPtsProf/NStdBins));
	  OutProfile.rstdq[i_bin] /= (b/(double)(NPtsProf/NStdBins));
	  OutProfile.rstdu[i_bin] /= (b/(double)(NPtsProf/NStdBins));
	  OutProfile.rstdv[i_bin] /= (b/(double)(NPtsProf/NStdBins));
	} 
      }
      else{	
	Zero2One(&OutProfile, NPtsProf, PSRName);
      }
      /* Scale profile to input template before outputting to file */
    }

  
    /* take average MJD of first to last scan */
    MJD_mid = (MJD_first + MJD_last)/2.;
    IMJDMid = floor(MJD_mid);
    MJDSecsMid = (MJD_mid - IMJDMid)*86400.;
  
    printf("MJD_mid = %lf, IMJDMid = %lf, MJDSecsMid = %lf\n",MJD_mid,IMJDMid,MJDSecsMid);fflush(stdout); 

    /* to choose a channel to put in the header for now, just use the average 
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
	    Hdr[0].obs.ObsvtyCode, 1, PSRName, 
	    0.);             
    //	  Subhdr.DumpRefPhase[i_chan]);             

    //}

  /* Open file for writing */
  if(Cmd->OutfileP)
    strcpy(Outfile, Cmd->Outfile);
  else
    sprintf(Outfile,"Template.out");
    
  /* now write the output ascii added profile */
  if ((Fout = fopen(Outfile,"w")) == 0)
    { printf("Cannot open %s. Exiting...\n",Outfile); exit(1); }

  fprintf(Fout,"%s\n",Header);

  if(Cmd->NoiseCutP) TemplateCutoff(&OutProfile, NStdBins, CutThresh);


  for(i_bin=0;i_bin<NPtsProf;i_bin++) {
    /* see how strong the linear polarization is */
    x = OutProfile.stdlin[i_bin]*OutProfile.Srms; 
    ptype = 43.1;
    if (x > 1.) ptype=43.2;
    if (x > 2.) ptype=43.3;
    if (x > 3.) ptype=43.4;
    if (x > 4.) ptype=43.5;
    if (x > 5.) ptype=43.6;
    fprintf(Fout,"%8d %22.7f %22.7f %22.7f %22.7f %22.7f %22.7f %22.7f %10.1f\n",i_bin,
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
  
  /* Output "SNR" and RMS of final profile */
  NStdBins = (int)NPtsProf; // New NStdBins
  Duty = DutyLookup(PSRName);
  BMask(OutProfile.rstds,&NStdBins,&Duty,FinalMask);
  Baseline(OutProfile.rstds,FinalMask,&NStdBins,&SBase,&Srms);
  SPeak =  FindPeak(OutProfile.rstds,&NStdBins,&spk);
  OutProfile.SNR = (SPeak-SBase)*Srms;

  printf("Output Template:  Peak = %.3lf,  Baseline RMS = %.3lf,   SNR = %.3lf\n\n",
	 SPeak, 1./Srms, OutProfile.SNR);
  
  printf("\nCompleted successfully.\n\n");fflush(stdout);
  
  exit(0);

}

/* Little subroutine that takes a profile and scales it to have a baseline of zero and a maximum value of 1 */
void Zero2One(struct StdProfs *Profile, int NBins, char *PSRName)
{

  int     i_bin, ngood, spk;
  float   profs[NBINMAX], amps[NBINMAX], phas[NBINMAX];
  float   Shift, EShift, SNR, ESNR, b, errb;
  double  SBase, QBase, UBase, VBase;
  double  Srms, Qrms, Urms, Vrms;
  double  Duty, SPeak, FinalMask[NBINMAX];

  cprofc(Profile->rstds, NBins, amps, phas);
  memcpy(profs, Profile->rstds, sizeof(float)*NBINMAX);
  /* Now fftfit to find scale factor */
  fftfit_(profs,&amps[1],&phas[1],
	  &NBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);
  /* Get Baseline value*/
  Duty = DutyLookup(PSRName);
  BMask(Profile->rstds, &NBins, &Duty, FinalMask);
  Baseline(Profile->rstds, FinalMask, &NBins, &SBase, &Srms);
  Baseline(Profile->rstdq, FinalMask, &NBins, &QBase, &Qrms);
  Baseline(Profile->rstdu, FinalMask, &NBins, &UBase, &Urms);
  Baseline(Profile->rstdv, FinalMask, &NBins, &VBase, &Vrms);
  SPeak =  FindPeak(Profile->rstds, &NBins, &spk); 

  //printf("%f,   %f,  %f,  %f,  %f\n", Duty, SBase, QBase, UBase, VBase);
  //printf("%f \n",amps[1]);

  /* Subtract baseline and scale by (peak-baseline) */
  /* Will remove baseline of each pol'n but scale factor is always b from 
    fftfit of total power profile */
  for(i_bin=0;i_bin<NBins;i_bin++) {
    Profile->rstds[i_bin] = (Profile->rstds[i_bin] - SBase)/(SPeak - SBase) ;
    Profile->rstdq[i_bin] = (Profile->rstdq[i_bin] - QBase)/(SPeak - SBase) ;
    Profile->rstdu[i_bin] = (Profile->rstdu[i_bin] - UBase)/(SPeak - SBase) ;
    Profile->rstdv[i_bin] = (Profile->rstdv[i_bin] - VBase)/(SPeak - SBase) ;
  }

}





void TemplateCutoff(struct StdProfs *StdProfile, int NBins, float CutThresh)
{

  int    i, j, k, kc, n2, noise_bin,  n_noise_bins, n_test;
  int    KCrit, KCritMin, KCritDiff; // critical harmonic bin number
  /*   int    forward=1,real=1; */
  int    back=-1,complex=0;
  float  MeanNoise=0;
  double f_W, f_B, FilterFunc[NBins/2], MinFunc;//, FilterFuncMed[NBins/2];
  double Diff, MaxDiff;
  FILE   *Famp, *Ffilt; //, *Ffiltmed;


  /* Fourier Transform the input profiles: */
  printf("Truncating noise in the Fourier domain...\n\n");

  /*   cprofc(StdProfile->rstds,NBins,
       StdProfile->stdamp,StdProfile->stdpha);  */
  cprofc(StdProfile->rstds,NBins,
	 StdProfile->stdamps,StdProfile->stdphas);
  cprofc(StdProfile->rstdq,NBins,
	 StdProfile->stdampq,StdProfile->stdphaq);
  cprofc(StdProfile->rstdu,NBins,
	 StdProfile->stdampu,StdProfile->stdphau);
  cprofc(StdProfile->rstdv,NBins,
	 StdProfile->stdampv,StdProfile->stdphav);


  /* Now do a minimization to figure out where the noise cutoff should be */

  /* First, estimate the noise power in each bin to be the mean power 
     level of the last one-quarter of the power spectrum */
  n_noise_bins = (int)floor(((double)NBins/2.) /4.);
  noise_bin = (NBins/2) - n_noise_bins;
  MeanNoise=0;
  n_test=0;
  for (k=noise_bin; k<NBins/2; k++){
    MeanNoise += StdProfile->stdamps[k];
    n_test++;
  }
  MeanNoise /= n_noise_bins;

  if((Ffilt = fopen("FilterFunc.dat", "w")) == NULL)
    { printf("Cannot open FilterFunc.dat. Exiting...\n"); exit(4); }   
  /*  if((Ffiltmed = fopen("FilterFuncMed.dat", "w")) == NULL)
      { printf("Cannot open FilterFuncMed.dat. Exiting...\n"); exit(4); }   */

  /* Now test critical harmonic kc across entire spectrum: */
  MinFunc = 1.0e12;
  MaxDiff = 0.0;
  KCritMin = 0;
  KCritDiff = 0;
  for(kc=0; kc<NBins/2; kc++){
    /* Now calculate Weber filter at each point in power spectum: */
    FilterFunc[kc] = 0.;
    for (k=0; k<NBins/2; k++){
      f_W = (double)( ( (StdProfile->stdamps[k])*(StdProfile->stdamps[k]) 
			- MeanNoise*MeanNoise ) /
		      ((StdProfile->stdamps[k])*(StdProfile->stdamps[k])) );
      /* Now assign 0 or 1 to "brick wall" filter */
      f_B = k>kc ? 0.0 : 1.0 ; 

      FilterFunc[kc] += (f_W - f_B)*(f_W - f_B);
    }
    if (FilterFunc[kc] < MinFunc) {
      KCritMin = kc;
      MinFunc = FilterFunc[kc];
    }
    Diff = fabs(FilterFunc[kc] - FilterFunc[kc-1]);
    if(kc > 0 && (Diff > MaxDiff)){
      MaxDiff = Diff;
      KCritDiff = kc;
    }

    fprintf(Ffilt, "%6d  %.5lf\n", kc, FilterFunc[kc]);
  }
  
  fclose(Ffilt);
  printf("KCritMin = %d\n", KCritMin);
  printf("KCritDiff = %d\n", KCritDiff);

  if (CutThresh < 0.0) {
    KCrit = KCritDiff;
    printf("Have chosen to truncate at the kc corresponding to the ");
    printf("largest difference between bins in minimizing function.\n");
  }
  else if (CutThresh == 0.0){
    KCrit = KCritMin;
    printf("Have chosen to truncate at the minimum k_c of the minimizing ");
    printf("function.\n");
  }
  else if (CutThresh == 1.0){
    KCrit = NBins/2;  // so that we'll never truncate anything
    printf("Have chosen not to truncate.\n");
  }
  else { // between 0.0 and 1.0
    KCrit = KCritMin + (int)(CutThresh * (float)((NBins/2) - KCritMin) );
    printf("Have chosen a cut threshold of %.2f, or %.1f%% of the way between ",
	   CutThresh, 100.*CutThresh);
    printf("the minimum k_c of the minimizing function and the maximum ");
    printf("k_c value\n");
  }
  
  printf("Threshold k_c = %d\n\n", KCrit);
  

  /* Now Median filter the Filter function to be minimized, 
     before minimization -- this will give a smoother function.  Here we 
     have used 5 nearest neighbour to calulate the median around each point  */
  //Median(FilterFunc, FilterFuncMed, NBins/2, 5);
  /* MinFunc = 1.0e12;
  for (kc=0; kc<NBins/2; kc++){
    if (FilterFuncMed[kc] < MinFunc) {
      KCritMed = kc;
      MinFunc = FilterFuncMed[kc];    
    }
    fprintf(Ffiltmed, "%6d  %.5lf\n", kc, FilterFuncMed[kc]);
    } */

  /* printf("noise_bin = %d, n_noise_bins = %d, n_test = %d\n", 
     noise_bin, n_noise_bins, n_test); */

  /* Output to file log(harmonic) vs. log(amplitude) for later plotting */

  if((Famp = fopen("ProfAmps.dat", "w")) == NULL) 
    { printf("Cannot open ProfAmps.dat. Exiting...\n"); exit(4); }   
  for(i=1;i<NBins/2;i++){ 
    //printf("%6d  %.5f   ", i, StdProfile->stdamps[i]);
    fprintf(Famp, "%6d  %.5f\n", i, StdProfile->stdamps[i]);
  }
  fclose(Famp);
  
  
  /* Now set all amplitudes in bins > kc equal to zero */
  for (k=KCrit; k<NBins/2; k++) StdProfile->stdamps[k] = 0.0;
  //  for (k=KCritDiff; k<NBins/2; k++) StdProfile->stdamps[k] = 0.0;
  //for (k=KCrit; k<NBins/2; k++) StdProfile->stdamps[k] = 0.0;

  /*    Pha1 = ToPhase;  */
    
  /*  for(i=1;i<NBins/2;i++) {
      StdProfile->stdphas[i] = fmod((StdProfile->stdphas[i]+(i)*Pha1),TWOPI);
      StdProfile->stdphaq[i] = fmod((StdProfile->stdphaq[i]+(i)*Pha1),TWOPI);
      StdProfile->stdphau[i] = fmod((StdProfile->stdphau[i]+(i)*Pha1),TWOPI);
      StdProfile->stdphav[i] = fmod((StdProfile->stdphav[i]+(i)*Pha1),TWOPI);

      } */


  /*   printf("totwgt: %f\n",*totwgt);   */
  uncprofc(StdProfile->stdamps,StdProfile->stdphas,NBins,
	   &StdProfile->stds[0][0]);
  uncprofc(StdProfile->stdampq,StdProfile->stdphaq,NBins,
	   &StdProfile->stdq[0][0]);
  uncprofc(StdProfile->stdampu,StdProfile->stdphau,NBins,
	   &StdProfile->stdu[0][0]);
  uncprofc(StdProfile->stdampv,StdProfile->stdphav,NBins,
	   &StdProfile->stdv[0][0]);

  ffft_(&StdProfile->stds[0][0], &NBins, &back, &complex);
  ffft_(&StdProfile->stdq[0][0], &NBins, &back, &complex);
  ffft_(&StdProfile->stdu[0][0], &NBins, &back, &complex);
  ffft_(&StdProfile->stdv[0][0], &NBins, &back, &complex);

  for(i=0;i<NBins;i++) 
    for(j=0;j<2;j++) {
      StdProfile->stds[i][j] /= (float)(NBins);
      StdProfile->stdq[i][j] /= (float)(NBins);
      StdProfile->stdu[i][j] /= (float)(NBins);
      StdProfile->stdv[i][j] /= (float)(NBins);
    }

  n2 = NBins*2;

  /*  ffft_(&StdProfile->stds[0][0], &NBins, &back, &complex);
      ffft_(&StdProfile->stdq[0][0], &NBins, &back, &complex);
      ffft_(&StdProfile->stdu[0][0], &NBins, &back, &complex);
      ffft_(&StdProfile->stdv[0][0], &NBins, &back, &complex); */

  for(i=0;i<NBins;i++) {
    StdProfile->rstds[i]= StdProfile->stds[i][0]/2.;
    StdProfile->rstdq[i]= StdProfile->stdq[i][0]/2.;
    StdProfile->rstdu[i]= StdProfile->stdu[i][0]/2.;
    StdProfile->rstdv[i]= StdProfile->stdv[i][0]/2.;
  }


  /**********************************************************************/

  printf("Noise truncation complete.\n\n");

}
