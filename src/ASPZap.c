/* Takes output stokes profile fits file and performs one of a choice of
 * diagnostics to zap appropriate scans/channels.  Zap information is 
 * placed into a "zap file" for inclusion on ASPFitsReader command
 * line for re-processing without zapped scans/channels in output file
 * 
 * -- R. Ferdman,  25 October 2010 */

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "fitsio.h"
#include "ZapCmdLine.h" 
#include "ASPCommon.h"
#include "polyco.h"
#include "cpgplot.h"
#include "cpgplot_def.h"


#define YES 1
#define NO  0

#define FREQMODE 0
#define TIMEMODE 1

#define YES_WRAP 1
#define NO_WRAP  0


/* Structure to bundle all grayscale properties into one object */
struct gray 
{
  int x_dim;        // Dimension of x-coordinate
  int y_dim;        // Dimension of y-coordinate
  int x_range[2];   // Min and max elements of x-coordinate to use
  int y_range[2];   // Min and max elements of y-coordinate to use
  float z[2];       // Min and max array values to use as background and 
                    // foreground, respectively
  float tr[6];      // Scaling array for axis labelling
  float *array;     // Array for plotting grayscale, to avoid overwriting original
};

/* Structure to bundle plotting preferences into on object */
struct plot_info 
{
  float x[2];         // Min and max x-axis values in plot
  float y[2];         // Min and max y-axis values in plot
  char x_label[64];   // x-axis label
  char y_label[64];   // y-axis label
};

// int MakePoly(Cmdline *, struct ASPHdr *);
int ReadZap(struct ASPHdr *, char *, int *);
int NormProf(int, float *, char *);
void GetGray(float *, int *, int, int, int, int, int, 
	     float, float, float, float,
	     struct gray *, struct plot_info *);
void PlotGray(struct gray *, float *, struct plot_info *);
void PlotProf(int, float *, float, float);
void GetHistArray(int, float *, int *, int *, float *);
int  GetHist(int, float *, int, float *, float *);
void PlotHist(int, float *, float *, char *, char *);
void PlotAllHist(int , float *, float *, float *, float *, float *, float *, 
		 float *, float *, int);
void GetHistLimits(float *, int);
void DrawHistLimits(int , float *, float *, float *);
void PlotAllProf(float *, int *, float *, int *, struct ASPHdr, float, float, float, float, 
		 struct gray *, struct plot_info *, struct gray *, struct plot_info *);
int  ProfZoom(float *, int *, struct ASPHdr, float *, int *, int);
int  ProfZap(float *, int *, struct ASPHdr, float *, int *, int *, int *, int);
int  UndoZap(float *, int *, struct ASPHdr, int *, int *, int *, int);
int  GetYesNo(char *);


int main(int argc, char *argv[])
{
  
  int first_pass=1, read_data=0;
  int NFirstTable, NumHDU, hdutype;
  struct ASPHdr Hdr;
  struct SubHdr SubHdr;

  struct StdProfs *Profile, TemplateProf;
  int OutRootIndex=0, LastSlashIndex=0, TotalRootIndex=0;
  char FitsFile[256], OutRoot[256], ProgName[32];
  char ZapFile[256], PazFile[256];
  double   **ASquared, **BSquared, **ReAconjB, **ImAconjB;
  long     **SampleCount;
  /* Dummy JyPerCounts for MakeStokes routine since not cal'ing */
  double   JyPerCount[4] = {1.0, 1.0, 1.0, 1.0};  
  fitsfile *Fin;
  FILE *FZap, *FPaz;

  int i, i_bin, i_chan, i_dump, i_array, status=0;
  int bin[NBINMAX], NBins, n_bin_hist;
  int fitsstatus=0;
  long NPtsProf=0;
  double RA, Dec;
  char NObs[3];
  char Headerline[256]; // header line in standard profile

  int   ngood, bad_array=0, bin_down=0, NHistArray; 
  int   *ProfWgt, *ProfChanMask, *ProfDumpMask, *TempChanMask, *TempDumpMask;
  float profs[NBINMAX],stdamps[NBINMAX], stdphas[NBINMAX], pha1; // for fftfit
  float Shift,EShift,SNR,ESNR,b,errb; // more fftfit variables
  double SBase,Srms;
  double Duty;
  double FinalMask[NBINMAX];
  double RMS;
  float *ProfRMS, *ProfShift, *ProfeShift, *ProfScale;
  float *ProfChan, *ProfDump, *ProfAll; //, *PhaseBin;;
  /* float ProfSNR; */
  float *HistArray;
  float *RMSBinVal, *RMSHist;
  float *ShiftBinVal, *ShiftHist, *eShiftBinVal, *eShiftHist, *ScaleBinVal, *ScaleHist;
  float Sideband, ChanBW;
  Cmdline *Cmd;

  double          StartMJD;
  struct Polyco   *Polycos;
  int             n_poly=0;

  int dev_rms_gray, dev_prof_gray, dev_prof_plot, dev_hist;
  int i_min, i_max;
  float xaxis_min, xaxis_max, yaxis_min, yaxis_max, x1, x2, y1, y2;
  struct gray RMSGray, ProfChanGray, ProfDumpGray;
  struct plot_info RMSInfo, ProfChanInfo, ProfDumpInfo;
  
  /* Interactive variables */
  int   quit=0, hist_quit=0, prof_quit=0, n_click=0;
  int   prof_mode=FREQMODE, reset_plot=1, final_check=0, cant_undo=1;
  int   first_zap=1, first_bad_dump=1, first_bad_chan=1;
  int   def_col = 0;
  int   *bad_dump, *bad_chan;
  int   ZoomIndex[2], ZapIndex[2], *ZapFlagChan, *ZapFlagDump;
  float x_input, y_input;
  float RMSLim[2], ShiftLim[2], eShiftLim[2], ScaleLim[2];
  float tempRMSLim[2], tempShiftLim[2], tempeShiftLim[2], tempScaleLim[2];
  float ClickVal[2];
  char  char_input[128], first_click_char[128];

  struct RunVars  RunMode;

                                                                             
 /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /* Store program name */
  strcpy(ProgName, argv[0]);

  strcpy(FitsFile, Cmd->Infile);
  
  LastSlashIndex = -1;
  for(i=0;i<strlen(FitsFile);i++){
    if(!strncmp(&FitsFile[i],"/",1))
      LastSlashIndex = i;
  }
  /* Get output file root name from input file name */
  OutRootIndex = -1;
  /*  for (i=strlen(FitsFile)-1;i>0;i--){
    if(!strncmp(&FitsFile[i],".",1)){
      OutRootIndex = i;
      break;
    }
    } */
  for (i=0;i<strlen(FitsFile);i++){
    if(!strncmp(&FitsFile[i],".asp",4) || !strncmp(&FitsFile[i],".fits",5)){
      OutRootIndex = i;
      break;
    }
  }
  if (OutRootIndex == -1) OutRootIndex = strlen(FitsFile) - 1;

  /*   for (i=0;i<strlen(FitsFile);i++){
    if(!(strcmp(&FitsFile[i],"stokes.fits"))){
      OutRootIndex = i-1;
      break;
    }
    } */

  if(OutRootIndex == 0){
    printf("WARNING: input file name does not follow .stokes.fits convention.\n");
    sprintf(OutRoot,"ASPZapOut");
  }
  else {
    TotalRootIndex = OutRootIndex-LastSlashIndex-1;
    strcpy(OutRoot,"\0");
    strncpy(OutRoot,&FitsFile[LastSlashIndex+1],TotalRootIndex);
    strcpy(&OutRoot[TotalRootIndex],"\0");
  }
  /* Find position of last "." */
  
  /* Open output zap file */
  sprintf(ZapFile, "%s.zap.dat", OutRoot);
  
  if(Cmd->ZapInP){
    if(!strcmp(ZapFile, Cmd->ZapIn))
      sprintf(ZapFile, "%s_new.zap.dat", OutRoot);
  }
  //  else {    
  if((FZap = fopen(ZapFile, "w")) == NULL){
    fprintf(stderr, "Error in opening output zap file %s.\n", ZapFile);
    exit(1);
  }
    //  }
  


  /* If paz-compatible output is requested, open data file for that. */
  if(Cmd->PazP){
    sprintf(PazFile, "%s.paz", OutRoot);
    if((FPaz = fopen(PazFile, "w")) == NULL){
      fprintf(stderr, "Error in opening output zap file %s.\n", PazFile);
      exit(1);
    }   
  }

  /* Open fits data file */

  if(fits_open_file(&Fin, FitsFile, READONLY, &status)){
    printf("Error opening FITS file %s !!!\n", FitsFile);
    exit(1);
  }
  
  /* Read in values for header variables */
  if(ReadHdr(&Hdr, Fin, &RunMode) < 0){
    printf("%s> Unable to read Header from file %s.  Exiting...\n",
	   ProgName,FitsFile);
    exit(2);
  }
  
  /* Dynamically allocate RunMode variables */
  if (AllocRunMode(&RunMode) < 0){
    printf("Could not allocate RunMode structure.  Exiting...\n");
    exit(2);
  }
  strcpy(RunMode.Infile,Cmd->Infile); 
  
  RunMode.NBins = RunMode.NBinsOut = Hdr.redn.RNBinTimeDump;
  

 printf("\n==========================\n");
  if(!strcmp(Hdr.gen.BEName, "xASP")) 
    printf("ASP FITS Header %s\n", Hdr.gen.HdrVer);
  else if(!strcmp(Hdr.gen.FitsType, "PSRFITS"))
    printf("PSRFITS Header %s\n", Hdr.gen.HdrVer);
  else
    printf("Unknown FITS Header %s\n", Hdr.gen.HdrVer);
  printf("==========================\n\n");fflush(stdout);
  
  printf("Input file:  %s\n\n",FitsFile);fflush(stdout);

  printf("PSR %s:\n",Hdr.target.PSRName);
  printf("--------------\n\n");
  printf("Centre Frequency: %6.1lf MHz\n\n",Hdr.obs.FSkyCent);fflush(stdout);
  //  if(Cmd->NoBaseP)
  //  printf("Baseline subtraction turned OFF.\n\n");fflush(stdout);

  /* If user gives number of histogram bins then use that, otherwise, 
     use hardcode for now */
  if(Cmd->NBinHistP)
    n_bin_hist = Cmd->NBinHist;
  else
    n_bin_hist = 25;

  /* Malloc histogram arrays */
  RMSBinVal = (float *)malloc(n_bin_hist*sizeof(float));
  RMSHist = (float *)malloc(n_bin_hist*sizeof(float));
  if(Cmd->TemplateP){
    /* Malloc histogram arrays */
    ShiftBinVal = (float *)malloc(n_bin_hist*sizeof(float));
    ShiftHist = (float *)malloc(n_bin_hist*sizeof(float));
    eShiftBinVal = (float *)malloc(n_bin_hist*sizeof(float));
    eShiftHist = (float *)malloc(n_bin_hist*sizeof(float));
    ScaleBinVal = (float *)malloc(n_bin_hist*sizeof(float));
    ScaleHist = (float *)malloc(n_bin_hist*sizeof(float));
  }



  // NDump = 0;
  /* Get number of HDUs in fits file */
  fits_get_num_hdus(Fin, &NumHDU, &status);
  RunMode.Dedisp = 0;
  /* Assume no reference frequency for dedispersion, i.e. infinite ref frequency */
  RunMode.DedispRefFreq = Hdr.obs.FSkyCent;
  /* Will always scale by off-pulse RMS since we are not using cals */
  RunMode.Scale = 1;
  RunMode.Swap = 0;
  RunMode.Verbose = Cmd->VerboseP;
  strcpy(RunMode.Source, Hdr.target.PSRName);


  /**** Get number of dumps based on which file type this is ****/  
  if(!strcmp(Hdr.gen.BEName, "xASP")) {

    RunMode.NDumps = Hdr.redn.RNTimeDumps;

#if 0
    if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
      NDump = NumHDU-3;  /* the "3" is temporary, depending on how 
			    many non-data tables we will be using */
    }
    else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
      NDump = (NumHDU-3)/2;
    }
    else{
      printf("Do not recognize FITS file version number in header.\n");
      printf("This header %s. Exiting...\n",Hdr.gen.HdrVer);fflush(stdout);
      exit(3);
    }

    /* Now check last dump to make sure it wrote properly if scan 
       was ctrl-c'd -- if so, reduce NDumps by 1 to avoid disaster */
    fits_get_num_hdus(Fin, &NumHDU, &status);
    fits_movabs_hdu(Fin,NumHDU,&hdutype, &status);
    /* find NPtsProf */
    fits_get_num_rows(Fin, &NPtsProf, &status);status=0;    
    if ((Hdr.redn.RNBinTimeDump!=(int)NPtsProf)){
      printf("Warning:  last dump did not write cleanly.\n");
      printf("Reducing number of dumps read in file by one...\n\n");
      fflush(stdout);
      NDump--;      
    }
    Hdr.redn.RNTimeDumps = NDump;
#endif

  }
  else {
    if(!strcmp(Hdr.gen.FitsType, "PSRFITS")) {
      /* Set to dedisperse input profiles before processing */
      if(!Cmd->NoDedispP)
		  RunMode.Dedisp = 1;
      else
	printf("Dedispersion turned off.\n\n");
      // NDump = Hdr.redn.RNTimeDumps;
      RunMode.NDumps = Hdr.redn.RNTimeDumps;
      NPtsProf = Hdr.redn.RNBinTimeDump;
    }
    else {
      /* Do not recognize data format! */
      fprintf(stderr, "ASPFitsReader ERROR: Unrecognized file format.\n");
      exit(1);
    }
  }
  //  if(Cmd->NFiles > 1)
  //    printf("Total number of dumps: %d\n",NDump);
  
  /* malloc pointers to input data arrays -- will use these to construct 
     Stokes parameters */
  ASquared    = (double **)malloc(Hdr.obs.NChan*sizeof(double));
  BSquared    = (double **)malloc(Hdr.obs.NChan*sizeof(double));
  ReAconjB    = (double **)malloc(Hdr.obs.NChan*sizeof(double));
  ImAconjB    = (double **)malloc(Hdr.obs.NChan*sizeof(double));
  SampleCount = (long    **)malloc(Hdr.obs.NChan*sizeof(long));
  
  /* for(i=0;i<Hdr.obs.NChan;i++){
    ASquared[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ASquared[i][0], NPtsProf);
    BSquared[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&BSquared[i][0], NPtsProf);
    ReAconjB[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ReAconjB[i][0], NPtsProf);
    ImAconjB[i]    = (double *)malloc(NPtsProf*sizeof(double));
    DZero(&ImAconjB[i][0], NPtsProf);
    SampleCount[i] = (long   *)malloc(NPtsProf*sizeof(long));
    LZero(&SampleCount[i][0], NPtsProf);
    } */


  /* malloc an array of size (Hdr.redn.RNTimeDumps * NChan) and initialize it 
     with ones to start */
  /* This will be our final "Zap" grid */
  ProfWgt = (int *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(int));
  for (i=0; i<Hdr.redn.RNTimeDumps*Hdr.obs.NChan; i++) ProfWgt[i]=1;
  

  /* If input zap file given, open it and read in to assign ProfWgt values
     accordingly */
  if(Cmd->ZapInP) {
    if(ReadZap(&Hdr, Cmd->ZapIn, ProfWgt) < 0){
      fprintf(stderr, "Problem reading input zap file or initializing zap ");
      fprintf(stderr, "array.  Exiting.");
      exit(2);
    }
  }
  


  printf("Number of channels:  %d\n",Hdr.obs.NChan);
  printf("Number of dumps:     %d\n\n",Hdr.redn.RNTimeDumps);

  Profile = (struct StdProfs *)malloc(Hdr.obs.NChan*sizeof(struct StdProfs));
  ProfAll = (float *)malloc(Hdr.redn.RNBinTimeDump*sizeof(float));
  // PhaseBin = (float *)malloc(Hdr.redn.RNBinTimeDump*sizeof(float));
  ProfChan = (float *)malloc(Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(float));
  ProfDump = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(float));
  ProfChanMask = (int *)malloc(Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
  ProfDumpMask = (int *)malloc(Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));


  /* Channel Bandwidth -- divide by NChan-1 since we are then in effect removing 
     half a channel at each end in the calculation, since the channel labels are 
     for the centres of the channels */
  ChanBW = fabsf((float)(Hdr.obs.ChanFreq[Hdr.obs.NChan-1] - Hdr.obs.ChanFreq[0]))/
    ((float)Hdr.obs.NChan -1);
  if(Hdr.obs.ChanFreq[Hdr.obs.NChan-1] > Hdr.obs.ChanFreq[0]) Sideband = 1.0;
  else Sideband = -1.0;

/* Conver Ra and Dec, mainly for use with parallactic angle calculation if needed */
  if(Hdr.target.RA > 0. && fabs(Hdr.target.Dec) > 0. ){
    RA = (Hdr.target.RA/24.)*TWOPI;
    Dec = Hdr.target.Dec*TWOPI/360.;
  }
  else{
    printf("\nRA and Dec are not in file.  Parallactic angle calculations\n");
    printf(" will be incorrect...\n");fflush(stdout);
  }
  /* Convert NObs to integer for possible use in calculating parallactic angle */
  sscanf(Hdr.obs.ObsvtyCode,"%s",NObs);


  /* Read/make polycos here, if requested by user, for eventual 
     realignment of profiles */
  
    /* Make polyco file on the fly if requested by user, which will be 
     seen and used by GetPoly routine */
  if (Cmd->ParFileP) {
    if(MakePoly(Cmd->ParFile, &Hdr) < 0){
      fprintf(stderr, "Could not make polycos. Exiting...\n");
      exit(2);
    }
    /* Finally, set appropriate variables to let main AFR program know that 
       it can go ahead and read in an existing par file (poly_final.dat in 
       this case) */
    Cmd->PolyfileP=1;
    Cmd->PolyfileC=1;
    Cmd->Polyfile = (char *) malloc(64);
    strcpy(Cmd->Polyfile,"poly_final.dat");
  }
  //printf("PolyfileP = %d, Polyfile = %s\n\n", Cmd->PolyfileP, Cmd->Polyfile);
 /* Read in polyco.dat file and get polyco structure if requested on 
     command line */
  if(Cmd->PolyfileP){
    if(Cmd->PolyfileC == 0) sprintf(Cmd->Polyfile,"polyco.dat");
    printf("Polyco file name:  %s\n",Cmd->Polyfile);
    StartMJD = (double)Hdr.obs.IMJDStart + 
      ((double)Hdr.obs.StartTime)/86400.;
    /* malloc Polyco structure to number of dumps */
    Polycos = (struct Polyco *)malloc(MAX_PC_SETS*Hdr.obs.NChan*
				      sizeof(struct Polyco));    
    for(i_chan=0;i_chan<Hdr.obs.NChan;i_chan++){
      if((n_poly=GetPoly(Cmd->Polyfile, RunMode.Source, 
			 &Polycos[i_chan*MAX_PC_SETS], 
			 Hdr.obs.ChanFreq[i_chan], 
			 StartMJD)) < 1) {
	printf("Could not find polycos for all input profiles of \n");
	printf("PSR %s given as input.  Will not shift profiles.\n",
	       RunMode.Source);
	Cmd->PolyfileP=0;
      }
      //printf("%lf MHz:  %d polyco sets found...\n", Hdr.obs.ChanFreq[i_chan], n_poly);
    }
    if (Cmd->PolyfileP) printf("Polycos successfully found.\n\n");
  }


  if(Cmd->ZapInP)
    printf("Will use input zap file %s to initialize zap array.\n\n", Cmd->ZapIn);



  /* printf("Rotating profile(s) by %6.3lf radians or %6.2lf degrees...\n\n",
     Cmd->ByAngle,Cmd->ByAngle*360./TWOPI); */

  printf("Writing omission data file %s\n\n",ZapFile);
  if(Cmd->PazP)  
    printf("Writing paz-compatible command line to file %s\n\n",PazFile);

  

  /* Malloc the temporary array to go into histogram making routine to be the maximum 
     possible, i.e. the same as each of the RMS, Shift, eShift, and Scale factor arrays */
  HistArray = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(float));

  /* Malloc the RMS array */
  ProfRMS = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(float));

  /*** Set up an array of shift/scale values wrt a template supplied by user, using fftfit ***/
  if(Cmd->TemplateP){
    
    /* Do this setup only the first time through */
    /* Read in template profile */
    if (ReadASPAsc(Cmd->Template, &Headerline[0], bin,  
		   &TemplateProf, &NBins) < 0) {
      printf("Error in reading file %s.\n",Cmd->Template);
      fflush(stdout);
      exit(1);
    }
    /* Check that template has same number of bins as the data file profiles */
    if(NBins != Hdr.redn.RNBinTimeDump){
      /* Bin down input profiles only for fftfitting by using a flag */
      bin_down = 1;
      printf("Template profile MUST have same number of bins as data profiles! (%d != %d)\n",
	     NBins, Hdr.redn.RNBinTimeDump);
      printf("Will bin down profiles for calculating cross-correlation shifts...\n");
      //      fflush(stdout);
      //      exit(1);
    }
    /* Now have standrad profile read in.  cprofc it: */
    cprofc(TemplateProf.rstds,NBins,TemplateProf.stdamps,TemplateProf.stdphas);
    
    /* Get amplitude and phase information from the standard profile */
    memcpy(stdamps,TemplateProf.stdamps,sizeof(float)*NBINMAX);
    memcpy(stdphas,TemplateProf.stdphas,sizeof(float)*NBINMAX);
    
    /* Rotate standard profile to zero phase unless user requests not to do so */
    pha1 = stdphas[1];
    for(i_bin=1;i_bin<(NBins/2)+1;i_bin++) 
      stdphas[i_bin] = fmod(stdphas[i_bin] -i_bin*pha1,TWOPI);
  

  /* Malloc the ProfShifts array */
    ProfShift = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(float));
    ProfeShift = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(float));
    ProfScale = (float *)malloc(Hdr.redn.RNTimeDumps*Hdr.obs.NChan*sizeof(float));
  }

  if(!strcmp(Hdr.gen.BEName, "xASP")) {
    /* Move to the first data table HDU if this is an ASP fits file */
    if(!strcmp(Hdr.gen.HdrVer,"Ver1.0"))
      fits_movnam_hdu(Fin, BINARY_TBL, "ASPOUT0", 0, &status);
    else if (!strcmp(Hdr.gen.HdrVer,"Ver1.0.1"))
      fits_movnam_hdu(Fin, ASCII_TBL, "DUMPREF0", 0, &status);
    
    /* Get the current HDU number */
    fits_get_hdu_num(Fin, &NFirstTable);
  }
  
  /* Only read data file if this is the first time through, or user wants a 
     re-read to see effect of zapping before quitting */
  while((first_pass || read_data) && !quit){
  
    FZero(ProfAll, Hdr.redn.RNBinTimeDump);
    FZero(ProfChan, Hdr.obs.NChan*Hdr.redn.RNBinTimeDump);
    IZero(ProfChanMask, Hdr.obs.NChan*Hdr.redn.RNBinTimeDump);
    FZero(ProfDump, Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump);
    IZero(ProfDumpMask, Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump);

  /* now start running through each dump */
  for(i_dump=0;i_dump<Hdr.redn.RNTimeDumps;i_dump++){


    /* If xASP format, move to next dump's data */
    if(!strcmp(Hdr.gen.BEName, "xASP")) {
      if(!strcmp(Hdr.gen.HdrVer,"Ver1.0")){
	fits_movabs_hdu(Fin, NFirstTable+i_dump, &hdutype, &status); 
      }
      else if(!strcmp(Hdr.gen.HdrVer,"Ver1.0.1")){
	fits_movabs_hdu(Fin,NFirstTable+(i_dump)*2+1,&hdutype,&status);
	fits_movrel_hdu(Fin, -1, NULL, &status);
      }
      fits_get_num_rows(Fin, &NPtsProf, &status);status=0; 
    }

    /* ReadASPStokes(&Hdr, &SubHdr, Fin, NPtsProf, 
       Profile, i_dump, Cmd->VerboseP); */

    /* Read in data arrays */
    if (ReadData(&Hdr, &SubHdr, &RunMode, Fin, i_dump,
		 NPtsProf, ASquared, BSquared, ReAconjB, ImAconjB, 
		 SampleCount) < 0){
      fprintf(stderr, "ASPZap ERROR: Could not read data from ");
      fprintf(stderr, "file %s (Error occured when attempting to read ",
	      Cmd->Infile);
      fprintf(stderr, "dump %d)", i_dump);
      exit(1);
    }
    

    for(i_chan=0;i_chan<Hdr.obs.NChan;i_chan++){
      
 	  /* Construct Stokes parameters */
	  MakeStokes(&Hdr, &RunMode, &Profile[i_chan], 
		     ASquared[i_chan], BSquared[i_chan],
		     ReAconjB[i_chan], ImAconjB[i_chan],
		     JyPerCount);
	  
	  /* Now dedisperse to the centre frequency before further 
	     processing, if required */	 
	  if(RunMode.Dedisp){
  	    if(i_dump==0 && i_chan==0) 
		  fprintf(stdout, "Dedispersing data...\n");
	    if (Dedisperse(&Profile[i_chan], &RunMode, 
			   &Hdr, &SubHdr, 
			   i_chan) < 0) {
	      fprintf(stderr,"ERROR: Cannot dedisperse data.\n");
	      exit(1);
	    }
	  }

	  /* Shift phases by appropriate amounts if requested on command line 
	     -- do for each dump and channel */
	  if(Cmd->PolyfileP){
	    if(i_dump==0 && i_chan==0) {
	      printf("Applying polyco-based phase shifts...\n");
	      fflush(stdout);
	    }
	    if(PhaseShift(&Polycos[i_chan*MAX_PC_SETS], n_poly, 
			  &Profile[i_chan], &RunMode,
			  &Hdr, &SubHdr, i_chan) < 0) {
	      printf("Unable to shift profile phases.  Exiting...\n");
	      fflush(stdout);
	      exit(11); 
	    }
	  }


	  i_array = i_chan*Hdr.redn.RNTimeDumps + i_dump;
	  
	  
	  /* Only calculate RMS and template cross correlation quantities is 
	     this is our first time through */
	  if(first_pass) {

	    /* Check that profile doesn't already have a zero weight attached 
	       to it (from and input zap file) already, so we don't waste time */
	    if(ProfWgt[i_array]!=0){
	    
	      /* Check that profile is not zeroed out, or contains NaNs */
	      bad_array = ArrayZero(Profile[i_chan].rstds, Hdr.redn.RNBinTimeDump);
	      
	      if (bad_array) {
		ProfWgt[i_array] = 0;
		if(bad_array==1) 
		  printf("Scan %d, channel %d is a zeroed-out profile\n", i_dump, i_chan);
		if(bad_array==2) 
		  printf("Scan %d, channel %d has at least some bins that are NaN\n",
			 i_dump, i_chan);		
		if(bad_array==3) 
		  printf("Scan %d, channel %d has at least some bins that are not finite\n",
			 i_dump, i_chan);		
	      }
	    }
	    // MakePol(&RunMode,RunMode.NBins, &Profile[i_chan]);
	    
	    Duty = DutyLookup(Hdr.target.PSRName);
	    BMask(Profile[i_chan].rstds,&Hdr.redn.RNBinTimeDump,&Duty,FinalMask);
	    Baseline(Profile[i_chan].rstds,FinalMask,&Hdr.redn.RNBinTimeDump,
		     &SBase,&Srms);
	    
	    
	    // SPeak =  FindPeak(Profile[i_chan].rstds,&Hdr.redn.RNBinTimeDump,&spk);
	    // Profile[i_chan].SNR = SPeak*Srms;
	    
	    /* Build up array of off-pulse RMSs -- 1D representation of a 2D array, for pgplot */
	    //i_array = i_dump*Hdr.obs.NChan + i_chan;
	    //  printf("i_array = %d, n_array = %d\n", i_array, Hdr.obs.NChan*Hdr.redn.RNTimeDumps);
	    //      if (Srms > 0.) 
	    
	    /* Can still have a bad profile that is not full of NaNs and is 
	       finite; if the signal is weak or if there is RFI such that the 
	       pulse is lost, BMask will have trouble creating a mask for 
	       getting off-pulse RMS.  FinalMask will be zeroes and RMS will be 
	       NaNs.  So in addition to testing for profiles flagged as being 
	       zeroed-out or having NaNs or not finite, we must test to see if 
	       the calculated RMS is reasonable, i.e. no itself a NaN and not 
	       infinite. */
	    RMS = (float)(1./Srms);
	    if(isnan(RMS) || !isfinite(RMS)) ProfWgt[i_array] = 0;
	    
	    if (ProfWgt[i_array])  // if the RMS is finite
	      ProfRMS[i_array] = RMS;
	    else 
	      ProfRMS[i_array] = 0.;  // will hopefully avoid NaNs this way...
	    
	    
	    /*      if(i_array==186){
		    printf("BLAH BLAH Srms[186] = %lf, ProfRMS[186] = %f\n", Srms, ProfRMS[186]);fflush(stdout);
		    double temp_count=0.0;
		    for(i=0;i<Hdr.redn.RNBinTimeDump;i++) {
		    
		    printf("%.3f %.3lf     ",Profile[i_chan].rstds[i], FinalMask[i]);
		    temp_count += FinalMask[i];
		    }
		    
		    printf("temp_count = %lf\nHOO-WAH!!!\n", temp_count)	 ; 
		    exit(0);
		    } */
	    
	    
	    //ProfRMS[i_array] = (float)(i_chan+i_dump);
	    //if(i_chan==8) ProfRMS[i_array] = 0.;
	    
	    /* If user has supplied a template profile, build up array of shifts, shift errors, 
	       and sclae factors */
	    if(Cmd->TemplateP){
	      struct StdProfs ShiftProf;
	      struct RunVars  ShiftRunMode;
	      
	      if(ProfWgt[i_array]) { // i.e. good data
		/* Bin down data profile here for shift calculation, if needed */
		if(bin_down){
		  ShiftRunMode.BinDown = 1;
		  ShiftRunMode.NBins = Hdr.redn.RNBinTimeDump;
		  ShiftRunMode.NBinsOut = NBins;
		  ShiftRunMode.Verbose = Cmd->VerboseP;
		  
		  BinDown(&ShiftRunMode, &Profile[i_chan], &ShiftProf);
		}
		else{
		  memcpy(&ShiftProf, &Profile[i_chan], 
			 sizeof(struct StdProfs));	  
		}
		memcpy(profs,ShiftProf.rstds,sizeof(float)*NBINMAX);
		/* Now fftfit to find shift required in second profile */
		fftfit_(profs,&stdamps[1],&stdphas[1],
			&NBins,&Shift,&EShift,&SNR,&ESNR,&b,&errb,&ngood);  
		ProfShift[i_array] = Shift;
		ProfeShift[i_array] = EShift;
		ProfScale[i_array] = b;
	      }
	    }
	    
	  }
	  /* Add prof in channels, and then in dumps, then all together, only using 
	     non-zapped scans */
	  if(ProfWgt[i_array]) { // i.e. good data
	    
	    for (i_bin=0; i_bin<Hdr.redn.RNBinTimeDump; i_bin++){
	      /* Total summed profile */
	      ProfAll[i_bin] += 
		((float)ProfWgt[i_array])*Profile[i_chan].rstds[i_bin];
	      /* profiles over frequencies  (NChan profiles) */
	      ProfChan[i_chan*Hdr.redn.RNBinTimeDump + i_bin] += 
		((float)ProfWgt[i_array])*Profile[i_chan].rstds[i_bin];
	      /* Profiles over time (NDump profiles) */
	      ProfDump[i_dump*Hdr.redn.RNBinTimeDump + i_bin] += 
		((float)ProfWgt[i_array])*Profile[i_chan].rstds[i_bin];
	      //  if(isnan(Profile[i_chan].rstds[i_bin]))
	      //  printf("***ALERT!!!*** i_dump = %d,  i_chan = %d, i_bin = %d is NaN!!\n", i_dump, i_chan, i_bin);
	      /* Make mask = 1 once the first non-zero-weight profile is found, since 
		 there will be at least *one* good profile going into the sum */
	      ProfChanMask[i_chan*Hdr.redn.RNBinTimeDump + i_bin] = 1;
	      ProfDumpMask[i_dump*Hdr.redn.RNBinTimeDump + i_bin] = 1;
	    }
	    /*	n_prof_all++;
		n_prof_chan[i_chan]++;
		n_prof_dump[i_dump]++; */
	  }
	  
	  /* Normalize each channel profile if this is the last dump */
	  if(i_dump == Hdr.redn.RNTimeDumps-1){
	    if(NormProf(Hdr.redn.RNBinTimeDump, 
			&ProfChan[i_chan*Hdr.redn.RNBinTimeDump], 
			&Hdr.target.PSRName[0]) < 0){
	      printf("Added Profile across channel %d (%.6lf) could not be ",
		     i_chan, Hdr.obs.ChanFreq[i_chan]);
	      printf("normalised.  Have left as is.\n");
	    }
	  }
	  
    }
    
    /* Normalize each dump profile now that it is all added at this point */
    if(NormProf(Hdr.redn.RNBinTimeDump, &ProfDump[i_dump*Hdr.redn.RNBinTimeDump], 
		Hdr.target.PSRName) < 0){
      printf("Added Profile across subint %d could not be ", i_dump);
      printf("normalised.  Have left as is.\n");
    }
    
    if(Cmd->VerboseP) printf("\n");fflush(stdout);
  }
  
  /* Normalize total added profile to have values between 0 and 1 */
  if (NormProf(Hdr.redn.RNBinTimeDump, ProfAll, Hdr.target.PSRName) < 0){
    printf("Total cumulative profile across all subints and channels could ");
    printf("not be normalised.  Have left as is.\n");
  }
  
  /* Finally, initialise current RMS, Shift, eShift, and Scale factor limits
   * before going interactive */
  if (first_pass) {
    tempRMSLim[0] = RMSLim[0] = FMin(ProfRMS, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_min);
    tempRMSLim[1] = RMSLim[1] = FMax(ProfRMS, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_max);
    if(Cmd->TemplateP){
      tempShiftLim[0] = ShiftLim[0] = FMin(ProfShift, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_min);
      tempShiftLim[1] = ShiftLim[1] = FMax(ProfShift, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_max);
      tempeShiftLim[0] = eShiftLim[0] = FMin(ProfeShift, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_min);
      tempeShiftLim[1] = eShiftLim[1] = FMax(ProfeShift, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_max);
      tempScaleLim[0] = ScaleLim[0] = FMin(ProfScale, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_min);
      tempScaleLim[1] = ScaleLim[1] = FMax(ProfScale, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_max);
    }
  }
  /* Done with FITS file.  Close it. */
  /*  fits_close_file(Fin, &fitsstatus); */



  if (first_pass) {
  /* Now set up pgplot for the three plot windows */
  /* RMS grayscale plot */
    if ((dev_rms_gray = cpgopen("/xs")) < 1) {
      printf("Warning:  pgplot device could not be initialized.\n");
      exit(1);
    }
//    cpgpap(8., 1.0);
    cpgpap(6.4, 1.0);
    cpgsch(1.);  
    
    /* Profile plots */
    if ((dev_prof_gray = cpgopen("/xs")) < 1) {
      printf("Warning:  pgplot device could not be initialized.\n");
      exit(1);
    }
    /* Make a long rectangular plot surface for profiles */
    cpgpap(3.3, 1.9);
//    cpgpap(5.5, 1.9);
	//  printf("VIEW WIDTH = %f,  VIEW HEIGHT = %f\n\n", view_width, view_height);
    cpgsch(1.);  
    
   if ((dev_prof_plot = cpgopen("/xs")) < 1) {
      printf("Warning:  pgplot device could not be initialized.\n");
      exit(1);
    }
    /* Make a long rectangular plot surface for profiles */
    cpgpap(5., 0.618034);
//    cpgpap(5., 0.618034);
    cpgsch(1.8);  


    /* Histogram plots */
    if ((dev_hist = cpgopen("/xs")) < 1) {
      printf("Warning:  pgplot device could not be initialized.\n");
      exit(1);
    }
    /* Make a long rectangular plot surface for 4 histograms if we are doing 
       template matching  */
    if(Cmd->TemplateP)
      cpgpap(2.6, 2.45);
//      cpgpap(4., 2.45);
    else
      cpgpap(4., 0.6180);
    cpgsch(1.9);  
    
  /* Set limiting x and y-axis coordinates for all the windows */
    xaxis_min = 0.11;
    xaxis_max = 0.97;
    yaxis_min = 0.12;
    yaxis_max = 0.97;

  }

  /* Now switch first_pass flag to note that we have gone through file once */
  first_pass = 0;
  /* Only redraw profiles if user asks (it's slow).  
     So assume this isn't the case */
  read_data = 0;

  /* continue until user quites */
  while (!read_data && !quit) {
    
  /**************** Plot RMS greyscale ******************/
  cpgslct(dev_rms_gray);


  cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max);
  sprintf(RMSInfo.x_label,"Integration number");
  sprintf(RMSInfo.y_label,"Observing frequency (MHz)");
  /* Add a half-channel at each end, to make labels refer to centred of channels */
  GetGray(ProfRMS, ProfWgt, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, 
	  1, Hdr.redn.RNTimeDumps, 1, Hdr.obs.NChan,
	  0., (float)(Hdr.redn.RNTimeDumps-1), 
	  Hdr.obs.ChanFreq[0], 
	  Hdr.obs.ChanFreq[Hdr.obs.NChan-1],
	  &RMSGray, &RMSInfo); 
  

  cpgeras();
  printf("Plotting RMS grayscale...\n");
  PlotGray(&RMSGray, ProfRMS, &RMSInfo);

  /******************************************************/


  /******************* Plot profiles *******************/
  cpgslct(dev_prof_gray);

  printf("Plotting Channel/Subint grayscale...\n");
  PlotAllProf(ProfChan, ProfChanMask, ProfDump, ProfDumpMask, Hdr,
	      xaxis_min, xaxis_max, yaxis_min, yaxis_max, 
	      &ProfChanGray, &ProfChanInfo, &ProfDumpGray, &ProfDumpInfo);
  
  cpgslct(dev_prof_plot);

  /* Now plot Added profile on the bottom */

  /* Bottom panel */
  x1 = xaxis_min;
  x2 = xaxis_max;
  y1 = yaxis_min;
  y2 = yaxis_max;

  printf("Plotting added pulse profile...\n");
  cpgsvp(x1, x2, y1, y2);

  PlotProf(Hdr.redn.RNBinTimeDump, ProfAll, 0.0, 1.0);

  /****************************************************/

  /**************** Plot histograms ******************/

  /* Now have all the arrays of diagnostic values we need.  Create histograms of each */
  GetHistArray(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfRMS, ProfWgt, &NHistArray, HistArray);
  //     printf("Made it IN!\n");fflush(stdout);
  //  if(GetHist(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfRMS, n_bin_hist, RMSBinVal, RMSHist) < 0){
  if(GetHist(NHistArray, HistArray, n_bin_hist, RMSBinVal, RMSHist) < 0){
    fprintf(stderr, "Error in forming histogram of RMS values.  Exiting.\n"); 
    exit(1);
  }
  //   printf("Made it OUT!\n");fflush(stdout);

  if(Cmd->TemplateP){
    
    GetHistArray(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfShift, ProfWgt, &NHistArray, HistArray);
    // if(GetHist(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfShift, n_bin_hist, ShiftBinVal, ShiftHist) < 0){
    if(GetHist(NHistArray, HistArray, n_bin_hist, ShiftBinVal, ShiftHist) < 0){
      fprintf(stderr, "Error in forming histogram of Shift values.  Exiting.\n"); 
      exit(1);
    }
    GetHistArray(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfeShift, ProfWgt, &NHistArray, HistArray);
    //    if(GetHist(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfeShift, n_bin_hist, eShiftBinVal, eShiftHist) <0){
    if(GetHist(NHistArray, HistArray,  n_bin_hist, eShiftBinVal, eShiftHist) <0){
      fprintf(stderr, "Error in forming histogram of eShift values.  Exiting.\n"); 
      exit(1);
    }
    GetHistArray(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfScale, ProfWgt, &NHistArray, HistArray);
    //    if(GetHist(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfScale, n_bin_hist, ScaleBinVal, ScaleHist) < 0){
    if(GetHist(NHistArray, HistArray, n_bin_hist, ScaleBinVal, ScaleHist) < 0){
      fprintf(stderr, "Error in forming histogram of Scale values.  Exiting.\n"); 
      exit(1);
    }
  }

  /* Plot entire set of histograms */
  cpgslct(dev_hist);
  cpgeras();
  /* Set limiting x and y-axis coordinates in the window */
  
 

 /* start with RMS plot, and if Template is given, do fftfit shift, shift errors, 
     and scaling factor */
  /* Set viewport coordinates for current plot */

  if(Cmd->TemplateP){
    //    cpgsvp(xaxis_min, xaxis_max, 0.03+yaxis_max-(yaxis_max-yaxis_min)/4., yaxis_max);
    cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max);
    cpgsubp(1,4);
  }
  else{
    cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max);
  }
  printf("Plotting histograms...\n");
  PlotAllHist(n_bin_hist, RMSBinVal, RMSHist, ShiftBinVal, ShiftHist, 
		   eShiftBinVal, eShiftHist, ScaleBinVal, ScaleHist, Cmd->TemplateP);
  

 /* Wait for user input */

  
  /* continue until user quites */
  // while (!quit) {

    /* if(!cpgcurs(&x_input, &y_input, char_input)){
      fprintf(stderr,"Error in accepting user input.\n");
      exit(1);     
      } 

    printf("x input = %f\n", x_input);
    printf("y input = %f\n", y_input);
    printf("char input = %s\n\n", char_input); */

    

    printf("\nEnter desired action: (h)istogram, (p)rofile, (r)e-read, or (q)uit:  ");
    fgets(&char_input[0], 64, stdin);
    
    /* Ensure that input is only one character in length */
    // if(srtrlem(char_input[0]) == 1){

      /* If user types 'h' that means they will want to adjust acceptable value 
	 limits based on histograms */
    //    if(!strncmp(char_input, "h\n", 2)) {
    if(!strncmp(char_input, "h", 1)) {
      hist_quit = 0;
      while (!hist_quit){
	/* Start with accepting input about min and max histogram ranges desired */
	cpgslct(dev_hist);
	//cpgsls(l_dashed); /* dashed line */
	cpgqci(&def_col); /* Save default colour, whatever that happens to be */
	cpgsci(c_red); /* make hashed regions red */
	
	
	if(Cmd->TemplateP) 
	  cpgpanl(1,1);
	printf("\nEnter desired limits for RMS histogram [leave as is]:  ");
	//printf("input string is _%s_ you know, and %d in length\n", input_str, (int)(strlen(input_str)));
	
	/* Get limits from user input */
	GetHistLimits(tempRMSLim, NO_WRAP);
	/* Draw lines to show chosen limits */
	DrawHistLimits(n_bin_hist, RMSBinVal, RMSHist, tempRMSLim);
	printf("You have chosen to exclude all profiles with off-pulse RMS outside ");
	printf("the limits [%.3f, %.3f]\n", 
	       tempRMSLim[0], tempRMSLim[1]);
	
	if(Cmd->TemplateP){
	  cpgpanl(1,2);
	  printf("\nEnter desired limits for cross correlation shift ");
	  printf("histogram [leave as is]:  ");
	  /* Get limits from user input */
	  GetHistLimits(tempShiftLim, YES_WRAP);
	  /* Draw lines to show chosen limits */
	  DrawHistLimits(n_bin_hist, ShiftBinVal, ShiftHist, tempShiftLim);
	  if(tempShiftLim[0] > tempShiftLim[1]){  /* i.e. will wrap around edge */
	    printf("You have chosen to exclude all profiles with shift *inside* ");
	    printf("the limits [%.1lf, %.1lf]\n", 
		   tempShiftLim[1], tempShiftLim[0]);
	  } 
	  else {
	    printf("You have chosen to exclude all profiles with shift *outside* ");
	    printf("the limits [%.1lf, %.1lf]\n", 
		   tempShiftLim[0], tempShiftLim[1]);
	  }
	  
	  cpgpanl(1,3);
	  printf("\nEnter desired x,y limits for cross correlation shift error  ");
	  printf("histogram [leave as is]:  ");
	  /* Get limits from user input */
	  GetHistLimits(tempeShiftLim, NO_WRAP);
	  /* Draw lines to show chosen limits */
	  DrawHistLimits(n_bin_hist, eShiftBinVal, eShiftHist, tempeShiftLim);
	  printf("You have chosen to exclude all profiles with shift error outside ");
	  printf("the limits [%.2lf, %.2lf]\n", 
		 tempeShiftLim[0], tempeShiftLim[1]);
	  
	  cpgpanl(1,4);
	  printf("\nEnter desired limits for cross correlation scale factor ");
	  printf("histogram [leave as is]:  ");
	  /* Get limits from user input */
	  GetHistLimits(tempScaleLim, NO_WRAP);
	  /* Draw lines to show chosen limits */
	  DrawHistLimits(n_bin_hist, ScaleBinVal, ScaleHist, tempScaleLim);
	  printf("You have chosen to exclude all profiles with scale factor ");
	  printf("outside the limits [%.2f, %.2f]\n", 
		 tempScaleLim[0], tempScaleLim[1]);
	  
	}
	cpgsci(def_col); /* restore default colour */
	
	if(GetYesNo("Are you satisfied with these histogram cuts?") == NO){
	  /* Redraw histogram plots and re-ask limits */
	  cpgslct(dev_hist);
	  PlotAllHist(n_bin_hist, RMSBinVal, RMSHist, ShiftBinVal, ShiftHist, 
		      eShiftBinVal, eShiftHist, ScaleBinVal, ScaleHist, Cmd->TemplateP);
	  /* Restore limits */
	  for (i=0; i<2; i++){
	    tempRMSLim[i] = RMSLim[i];
	    if(Cmd->TemplateP){
	      tempShiftLim[i] = ShiftLim[i];
	      tempeShiftLim[i] = eShiftLim[i];
	      tempScaleLim[i] = ScaleLim[i];
	    }
	  }
	  
	}
	else {  /* Must be yes since can't get out of GetYesNo otherwise */
	  for (i=0; i<2; i++){
	    /* Make histogram cuts permanent */
	    RMSLim[i] = tempRMSLim[i];
	    if(Cmd->TemplateP){
	      ShiftLim[i] = tempShiftLim[i];
	      eShiftLim[i] = tempeShiftLim[i];
	      ScaleLim[i] = tempScaleLim[i];
	    }
	  }
	  /* exit input of histograms, set various Wgt array elements to zero, and
	     re-plot everything for next round. */
	  hist_quit = 1;
	}
      }
      printf("Will now adjust plots to reflect requested zapping...\n");
      
      for(i_array=0;i_array<Hdr.redn.RNTimeDumps*Hdr.obs.NChan;i_array++){
	if(ProfRMS[i_array] < RMSLim[0] || 
	   ProfRMS[i_array] > RMSLim[1]) 
	  ProfWgt[i_array] = 0;
	if(Cmd->TemplateP){
	  if(ShiftLim[0] < ShiftLim[1]){
	    if(ProfShift[i_array] < ShiftLim[0] || 
	       ProfShift[i_array] > ShiftLim[1]) 
	      ProfWgt[i_array] = 0;
	  }
	  else{
	    if(ProfShift[i_array] < ShiftLim[0] && 
	       ProfShift[i_array] > ShiftLim[1]) 
	      ProfWgt[i_array] = 0;
	  }
	  if(ProfeShift[i_array] < eShiftLim[0] || 
	     ProfeShift[i_array] > eShiftLim[1]) 
	    ProfWgt[i_array] = 0;
	  if(ProfScale[i_array] < ScaleLim[0] || 
	     ProfScale[i_array] > ScaleLim[1]) 
	    ProfWgt[i_array] = 0;
	}
	ProfRMS[i_array] = (float)ProfWgt[i_array] * ProfRMS[i_array];
      }
      
      /* Return colours and line styles back to how they were before */
      //cpgsls(l_full); /* normal solid line */
      
    }
    //    else if(!strncmp(char_input, "p\n", 2)) {
    else if(!strncmp(char_input, "p", 1)) {

      printf("Will now perform cuts based on the frequency and/or ");
      printf("time-added profiles...\n\n");

      /* Firstly, make a temporary Mask for each of Chan and Dump arrays */
      TempChanMask = (int *)malloc(Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
      TempDumpMask = (int *)malloc(Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));
      memcpy(&TempChanMask[0], &ProfChanMask[0], 
	     Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
      memcpy(&TempDumpMask[0], &ProfDumpMask[0], 
	     Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));

      /* Set up ZapFlags, and set contents to zero */
      ZapFlagChan = (int *)malloc(Hdr.obs.NChan*sizeof(int));
      ZapFlagDump = (int *)malloc(Hdr.redn.RNTimeDumps*sizeof(int));
      IZero(ZapFlagChan, Hdr.obs.NChan);
      IZero(ZapFlagDump, Hdr.redn.RNTimeDumps);

      cpgslct(dev_prof_gray);
      /*      cpgeras();
      sprintf(ProfChanInfo.x_label,"Profile phase");
      sprintf(ProfChanInfo.y_label,"Observing frequency (MHz)");
      cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max); 
      PlotGray(ProfChanGray, ProfChan, ProfChanInfo); */
      /* User input now will be mouse clicks on the profile pgplot window */

      /***      if(!cpgband(5, 0, 0., 0., &x_input, &y_input, char_input)){
		fprintf(stderr,"Error in accepting user input.\n");
		exit(1);     
	}  ***/

      /*      printf("x input = %f\n", x_input);
	      printf("y input = %f\n", y_input);
	      printf("char input = %s\n\n", char_input);  */
      
      /* Start by plotting frequency mode only */
      prof_mode = FREQMODE;
      prof_quit = 0;
      n_click = 0;
      reset_plot = 1;
      final_check = 0;
      cant_undo = 1;
      /* Run a while loop for user input as with histogram choice */
      while (!prof_quit){
	
	/* First check whether previous command involves resetting profile */
	if(reset_plot){
	  cpgeras();
	  if(prof_mode==FREQMODE){
	    /* Re-plot profile window to be only frequency-added */
	    GetGray(ProfChan, TempChanMask, Hdr.redn.RNBinTimeDump*Hdr.obs.NChan, 
		    1, Hdr.redn.RNBinTimeDump, 1, Hdr.obs.NChan, 
		    0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
		    Hdr.obs.ChanFreq[0], // - Sideband*0.5*ChanBW, 
		    Hdr.obs.ChanFreq[Hdr.obs.NChan-1], // + Sideband*0.5*ChanBW,
		    &ProfChanGray, &ProfChanInfo);
	    
	    cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max); 
	    PlotGray(&ProfChanGray, ProfChan, &ProfChanInfo);
	    ZoomIndex[0] = 0;
	    ZoomIndex[1] = Hdr.obs.NChan - 1;
	  }
	  else if (prof_mode==TIMEMODE){
	    /* Re-plot profile window to be only frequency-added */
	    GetGray(ProfDump, TempDumpMask, Hdr.redn.RNBinTimeDump*Hdr.redn.RNTimeDumps, 
		    1, Hdr.redn.RNBinTimeDump, 1, Hdr.redn.RNTimeDumps, 
		    0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
		    0.0, (float)(Hdr.redn.RNTimeDumps-1),
		    &ProfDumpGray, &ProfDumpInfo);
	    	    
	    cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max); 
	    PlotGray(&ProfDumpGray, ProfDump, &ProfDumpInfo);
	    ZoomIndex[0] = 0;
	    ZoomIndex[1] = Hdr.redn.RNTimeDumps - 1;
	  }
	  n_click = 0;	  
	  reset_plot = 0;
	}

	if(final_check){
	  final_check = 0;
	  /* Plot both grayscales side by side again with current zaps */
	  PlotAllProf(ProfChan, TempChanMask, ProfDump, TempDumpMask, Hdr,
		      xaxis_min, xaxis_max, yaxis_min, yaxis_max, 
		      &ProfChanGray, &ProfChanInfo, &ProfDumpGray, &ProfDumpInfo);
	  if(GetYesNo("Are you satisfied with these frequency/subint cuts?") == YES){
	    printf("Will now make these cuts permanent and add them to current zap list.\n");
	    /* Now copy TempMask into the official mask */
	    memcpy(&ProfChanMask[0], &TempChanMask[0], 
		   Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
	    memcpy(&ProfDumpMask[0], &TempDumpMask[0], 
		   Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));

	    /* Now copy indices in channel and dump arrays to reflect in the ProfWgt array */
	    for(i_dump=0; i_dump<Hdr.redn.RNTimeDumps; i_dump++){
	      for(i_chan=0; i_chan<Hdr.obs.NChan; i_chan++){
		i_array = i_chan*Hdr.redn.RNTimeDumps + i_dump;
		if(ZapFlagDump[i_dump] || ZapFlagChan[i_chan]) ProfWgt[i_array]=0;
	      }
	    }


	    prof_quit = 1;
	  }
	  else{
	    printf("Current round of cuts have been reset. Restarting profile cuts.\n");
	    memcpy(&TempChanMask[0], &ProfChanMask[0], 
		   Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
	    memcpy(&TempDumpMask[0], &ProfDumpMask[0], 
		   Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));
	    prof_mode = FREQMODE;
	    reset_plot = 1;
	    cant_undo = 1;
	  }
	}
	else{

	/* Now ask for user input on the profile pgplot window */
	if(!cpgband(5, 0, 0., 0., &x_input, &y_input, char_input)){
	  fprintf(stderr,"Error in accepting user input.\n");
	  exit(1);     
	} 
	n_click++;



	/* Now start seeing what input values are */
	/* First choice is to see if user quits out of profile zapping ('q') */
	//	if (!strncmp(char_input, "q\n", 2)){
	if (!strncmp(char_input, "q", 1)){
	    final_check = 1;
	    reset_plot = 1;
	}
	/* Profile mode help ('?') */
	//	if (!strncmp(char_input, "?\n", 2)){
	if (!strncmp(char_input, "?", 1)){
	  printf("-------------------------------------------------------------------------\n");
	  printf("Profile mode command options:\n");
	  printf("-------------------------------------------------------------------------\n");
	  printf("    Left click - left click:  Zoom into region bounded by clicks\n");
	  printf("                          u:  Undo zoom; regain full view range of \n");
	  printf("                              channels/subints\n");
	  printf("   Left click - right click:  Zap range of channels/subints bounded by ");
	  printf("                              mouse clicks\n");
	  printf("  Right click - right click:  Undo previous zap\n");
	  printf("                          q:  Quit profile mode; will be prompted to");
	  printf("                              verify zap choices made.\n");
	  printf("                          ?:  This help\n");
	  printf("-------------------------------------------------------------------------\n\n");
	  n_click = 0;
	}
	/* Let user pick frequency or time-added profile ('f' or 't').  If so, reset 
	   corresponding plot */
	//	else if(!strncmp(char_input, "f\n", 2)){
	else if(!strncmp(char_input, "f", 1)){
	  prof_mode = FREQMODE;
	  reset_plot = 1;
	  cant_undo = 1;
	}
	//	else if(!strncmp(char_input, "t\n", 2)){
	else if(!strncmp(char_input, "t", 1)){
	  prof_mode = TIMEMODE;
	  reset_plot = 1;
	  cant_undo = 1;
	}
	/* Now allow user to undo zoom with 'u'.  Simple flag switch, relying on 
	   prof_mode value to reset to correct plot */
	//	else if(!strncmp(char_input, "u\n", 2)){
	else if(!strncmp(char_input, "u", 1)){
	  reset_plot = 1;
	}

	/* Now start clicking options.   */
	/* For both cases, have one action if this is firt click, 
	   and another if it is the second click. */

	/* Left-click */
	//	else if(!strncmp(char_input, "A\n", 2)){
	else if(!strncmp(char_input, "A", 1)){
	  /* Set first zoom range value to this click */
	  if(n_click==1){
	    strcpy(first_click_char, char_input);
	    ClickVal[0] = y_input;
	  }
	  else{
	    /* Now we are on the second click. First set click value into 
	       increasing numerical value */
	    if(y_input > ClickVal[0]){
	      ClickVal[1] = y_input;
	    }
	    else{
	      ClickVal[1] = ClickVal[0];
	      ClickVal[0] = y_input;
	    }

	    /* If first click was left-click, then this means to zoom.  Otherwise, 
	       it is an unrecognized command, and we reset clicks, leaving 
	       current plot untouched. */
	    //	    if(!strncmp(first_click_char, "A\n", 2)) {
	    if(!strncmp(first_click_char, "A", 1)) {
	      if(prof_mode==FREQMODE){
		if(ProfZoom(ProfChan, TempChanMask, Hdr, ClickVal, ZoomIndex, prof_mode) < 0){
		  fprintf(stderr, "Could not zoom into desired frequency limits. Exiting.\n");
		  exit(2);
		}
	      }
	      else{
		if(ProfZoom(ProfDump, TempDumpMask, Hdr, ClickVal, ZoomIndex, prof_mode) < 0){
		  fprintf(stderr, "Could not zoom into desired subint limits. Exiting.\n");
		  exit(2);
		}
	      }
	      
	    }
	    else{
	      printf("Unrecognized command %s. Try again.\n", first_click_char);
	    }
	    n_click=0;
	  }
	
	  

	}

	/* Right-click */
	//	else if(!strncmp(char_input, "X\n", 2)){
	else if(!strncmp(char_input, "X", 1)){
	  /* Set first zoom range value to this click */
	  if(n_click==1){
	    strcpy(first_click_char, char_input);
	    ClickVal[0] = y_input;
	  }
	  else{
	    /* Now we are on the second click. First set click value into 
	       increasing numerical value */
	    if(y_input > ClickVal[0]){
	      ClickVal[1] = y_input;
	    }
	    else{
	      ClickVal[1] = ClickVal[0];
	      ClickVal[0] = y_input;
	    }

	    /* If first click was left-click, then this means to zap. */
	    //	    if(!strncmp(first_click_char, "A\n", 2)) {
	    if(!strncmp(first_click_char, "A", 1)) {
	      
	      if(prof_mode==FREQMODE){
		if(ProfZap(ProfChan, TempChanMask, Hdr, ClickVal, ZoomIndex, 
			   ZapIndex, ZapFlagChan,  prof_mode) < 0){
		  fprintf(stderr, "Could not zap desired frequency range. Exiting.\n");
		  exit(2);
		}	
	      }
	      else{
		if(ProfZap(ProfDump, TempDumpMask, Hdr, ClickVal, ZoomIndex, 
			   ZapIndex, ZapFlagDump, prof_mode) < 0){
		  fprintf(stderr, "Could not zap desired subint range. Exiting.\n");
		  exit(2);		
		}

	      }
	      cant_undo = 0;
	    }
	    /* If first click was right click, then this means to un-zap last zap range. */
	    //	    else if (!strncmp(first_click_char, "X\n", 2)) {
	    else if (!strncmp(first_click_char, "X", 1)) {
	      if(cant_undo){
		printf("Cannot undo last zap.  Either there is nothing to undo, ");
		printf("or you have already undone previous zap, or you have switched ");
		printf("from frequency to subint plot.\n");
	      }
	      else {
		if(prof_mode==FREQMODE){
		  if(UndoZap(ProfChan, TempChanMask, Hdr, ZoomIndex, ZapIndex, ZapFlagChan, 
			     prof_mode) < 0){
		    fprintf(stderr, "Could not undo last zap. Exiting.\n");
		    exit(2);		
		  }
		}
		else{
		  if(UndoZap(ProfDump, TempDumpMask, Hdr, ZoomIndex, ZapIndex, ZapFlagDump, 
			     prof_mode) < 0){
		    fprintf(stderr, "Could not undo last zap. Exiting.\n");
		    exit(2);		
		  }
		}
		printf("You have undone last zap\n");
		cant_undo=1;
	      }
	    }
	    /* Otherwise, it is an unrecognized command, and we reset clicks, leaving 
	       current plot untouched.  */
	    else{
	      printf("Unrecognized command %s. Try again.\n", first_click_char);
	    }
	    n_click=0;
	  }

	}

	/* If user has not entered one of the available options, then just reset click number, 
	   leaving current plot untouched. */
	else{
	  printf("Unrecognized command %s. Try again.\n", char_input);
	  n_click=0;
	}


	/*	printf("x input = %f\n", x_input);
		printf("y input = %f\n", y_input);
		printf("char input = %s\n\n", char_input);  */
			
	}
      /* Take note of zapped frequencies/subints */
      }
      
      free(TempChanMask);
      free(TempDumpMask);
      free(ZapFlagChan);
      free(ZapFlagDump);

    }
    //    else if(!strncmp(char_input, "r\n", 2)) {
    else if(!strncmp(char_input, "r", 1)) {
      read_data=1;
      printf("Re-reading data file, updating zap mask...\n\n");
    }
    //    else if(!strncmp(char_input, "q\n", 2)) {
    else if(!strncmp(char_input, "q", 1)) {
      quit=1;
    }
    else{
      printf("Unrecognised input.\n\n");
    }
    
    
  }
  }
  
  /* First look at ProfWgt and determine if there is an entire scan or channel 
     has been zapped */
  
  bad_dump = (int *)malloc(Hdr.redn.RNTimeDumps*sizeof(int));
  bad_chan = (int *)malloc(Hdr.obs.NChan*sizeof(int));
  IZero(bad_dump, Hdr.redn.RNTimeDumps);
  IZero(bad_chan, Hdr.obs.NChan);
  
  for (i_dump=0; i_dump<Hdr.redn.RNTimeDumps; i_dump++){
    for (i_chan=0; i_chan<Hdr.obs.NChan; i_chan++){
      i_array = i_chan*Hdr.redn.RNTimeDumps + i_dump;
      if(ProfWgt[i_array]==0){
	bad_dump[i_dump]++;
	bad_chan[i_chan]++;
      }
    }
  }
  
  /* Now go through and for each entire dump or frequency channel identified as bad, 
     write line out to file */

  for (i_dump=0; i_dump<Hdr.redn.RNTimeDumps; i_dump++){
    if(bad_dump[i_dump]==Hdr.obs.NChan){
      bad_dump[i_dump]=1;
      /* Now write this to file */
      fprintf(FZap, "  %d  -1.0\n", i_dump);
      if(Cmd->PazP){
	if(first_zap) {
	  fprintf(FPaz, "paz ");
	  first_zap = 0;
	}
	if(first_bad_dump) {
	  fprintf(FPaz, "-w \"");
	  first_bad_dump = 0;
	}
	fprintf(FPaz, "%d ", i_dump);
      }
    }
    else{
      bad_dump[i_dump]=0;
    }
  }
  if(Cmd->PazP && !first_bad_dump) fprintf(FPaz, "\" "); /* Close off last bad dump */

  for (i_chan=0; i_chan<Hdr.obs.NChan; i_chan++){
    if(bad_chan[i_chan]==Hdr.redn.RNTimeDumps){
      bad_chan[i_chan] = 1;
      /* Now write this to file */
      fprintf(FZap, "  -1  %.8lf\n", Hdr.obs.ChanFreq[i_chan]);
      if(Cmd->PazP){
	if(first_zap) {
	  fprintf(FPaz, "paz ");
	  first_zap = 0;
	}
	if(first_bad_chan) {
	  fprintf(FPaz, "-z \"");
	  first_bad_chan = 0;
	}
	fprintf(FPaz, "%d ", i_chan);
      }      
    }
    else{
      bad_chan[i_chan] = 0;
    }
  }
  if(Cmd->PazP && !first_bad_chan) fprintf(FPaz, "\" "); /* Close off last bad chan */
 
  /* Now go through each ProfWgt=0 scan/chan combination and if it is not part of
     an identified entire bad scan or channel, write it to file */ 
  for (i_dump=0; i_dump<Hdr.redn.RNTimeDumps; i_dump++){
    if(!bad_dump[i_dump]){
      for (i_chan=0; i_chan<Hdr.obs.NChan; i_chan++){
	if(!bad_chan[i_chan]){
	  /* Write individual scan/chan combo to file */
	  i_array = i_chan*Hdr.redn.RNTimeDumps + i_dump;
	  if (ProfWgt[i_array]==0)
	    fprintf(FZap, "  %d  %.8lf\n", i_dump, Hdr.obs.ChanFreq[i_chan]);	  
	  if(Cmd->PazP){
	    if(first_zap) {
	      fprintf(FPaz, "paz ");
	      first_zap = 0;
	    }
	    fprintf(FPaz, "-I %d %d ", i_dump, i_chan);
    	  }
	}
	
      }
    }
  }

  fclose(FZap);

  if(Cmd->PazP){
    fprintf(FPaz, "%s", FitsFile);
    fclose(FPaz);
  }
  
  /* Free ASquared, BSquared, etc. */
  free(ASquared);
  free(BSquared);
  free(ReAconjB);
  free(ImAconjB);
  free(SampleCount);

  free(RMSBinVal);
  free(RMSHist);
  if(Cmd->TemplateP){
    free(ShiftBinVal);
    free(ShiftHist);
    free(eShiftBinVal);
    free(eShiftHist);
    free(ScaleBinVal);
    free(ScaleHist);
  }
  /******************************************************/
  
  /* Done with FITS file.  Close it. */
  fits_close_file(Fin, &fitsstatus);

  
  printf("Completed successfully.\n\n");fflush(stdout);
  
  exit(0);
  
}

/* Normalize profile to have values between zero and one */
int NormProf(int NBin, float *Profile, char *Source)
{
 
  int i_bin, i_spk;
  double Duty,SPeak;
  double FinalMask[NBINMAX];
  double SBase, Srms; 
  
  /* Take baseline */
  Duty = DutyLookup(Source);
  BMask(Profile, &NBin, &Duty, FinalMask);
  Baseline(Profile, FinalMask, &NBin, &SBase,&Srms);
  SPeak =  FindPeak(Profile, &NBin, &i_spk);
  
  /* Safeguard against division by zero */
  if(SPeak == SBase) return -1;

  //  ProfSNR = (float)(SPeak*Srms);
  /* Normalise to be between 0. and 1. */
  for(i_bin=0;i_bin<NBin;i_bin++) 
    Profile[i_bin] = (Profile[i_bin] - (float)SBase)/(float)(SPeak - SBase);

  //  printf("HELLS YEAH:  SBase = %lf,  SPeak = %lf\n", SBase, SPeak);
  return 1;

}

/* "Nice and neat" routine to package all the plotting info for grayscale plots */
void GetGray(float *array, int *mask, int n_array, 
	     int i_xmin, int i_xmax, 
	     int i_ymin, int i_ymax, 
	     float x_min, float x_max, float y_min, float y_max,
	     struct gray *Gray, struct plot_info *PlotInfo)
{
  /* Covers all the quantities that would change between incarnations of this plot */
  int i_array, n_mask=0, i_min, i_max;
  float *mask_array;

  Gray->array = (float *)malloc(n_array*sizeof(float));
  memcpy(Gray->array, array, n_array*sizeof(float));
  mask_array = (float *)malloc(n_array*sizeof(float));

  Gray->x_dim = i_xmax - i_xmin + 1;
  Gray->y_dim = i_ymax - i_ymin + 1;

  Gray->tr[1] = (x_max - x_min)/(float)(Gray->x_dim-1);
  Gray->tr[2] = 0.0;
  Gray->tr[0] = x_min - Gray->tr[1];
  Gray->tr[4] = 0.0;
  Gray->tr[5] = (y_max - y_min)/(float)(Gray->y_dim-1);
  Gray->tr[3] = y_min - Gray->tr[5];

  /* Get z values from only non-masked parts of main array */
  for (i_array=0; i_array<n_array; i_array++)
    if(mask[i_array]) mask_array[n_mask++] = array[i_array];
  
  //  Gray->z[1] = FMax(array, n_array, &i_max);
  //  Gray->z[0] = FMin(array, n_array, &i_min);
  Gray->z[1] = FMax(mask_array, n_mask, &i_max);
  Gray->z[0] = FMin(mask_array, n_mask, &i_min);

  /* Make bad array values equal to less than the limit in z-values, 
     so it will show up black */
  for (i_array=0; i_array<n_array; i_array++)
    if(!mask[i_array]) Gray->array[i_array] = Gray->z[0] - 1.0;

  Gray->x_range[0] = i_xmin;
  Gray->x_range[1] = i_xmax;
  Gray->y_range[0] = i_ymin;
  Gray->y_range[1] = i_ymax;
  
  PlotInfo->x[0] = x_min - 0.5*Gray->tr[1];
  PlotInfo->x[1] = x_max + 0.5*Gray->tr[1];
  PlotInfo->y[0] = y_min - 0.5*Gray->tr[5];
  PlotInfo->y[1] = y_max + 0.5*Gray->tr[5];

  free(mask_array);

}



void PlotGray(struct gray *Gray, float *array, struct plot_info *PlotInfo) 
{

  cpgbbuf();  
  /* Set world coordinates for current plot */
  cpgswin(PlotInfo->x[0], PlotInfo->x[1], PlotInfo->y[0], PlotInfo->y[1]);
  /* Do grayscale */
  cpggray(Gray->array, Gray->x_dim, Gray->y_dim, Gray->x_range[0], Gray->x_range[1], 
	  Gray->y_range[0], Gray->y_range[1], Gray->z[1], Gray->z[0], Gray->tr);

  /* Set up plot and labelling preferences */
  cpgbox("BCTSN", 0.0, 0.0, "BCTSVN", 0.0, 0.0);            
  cpgmtxt("B", 2.35, 0.5, 0.5, PlotInfo->x_label);
  cpgmtxt("L", 3.1, 0.5, 0.5, PlotInfo->y_label);
  cpgebuf();  

}

/* Plot 1-D pulse profile within phase range [x1, x2] */
void PlotProf(int NBin, float *ProfArray, float x1, float x2) 
{
  
  int   i_bin, min_bin, max_bin, NewBins, i_min, i_max;
  // static int i_plot=0;
  float *PhaseBin, *Profile;
  char  x_label[64], y_label[64];

  if(x1 > x2){
    fprintf(stderr, "Choice of minimum profile phase (%.2f) is greater than or ", x1);
    fprintf(stderr, "equal to maximum profile phase (%.2f).\n", x2);
    exit(1);
  }

  min_bin = (int)(x1*((float)(NBin-1)));
  max_bin = (int)(x2*((float)(NBin-1)));

  NewBins = max_bin-min_bin+1;

  PhaseBin = (float *)malloc(NewBins*sizeof(float));
  Profile = (float *)malloc(NewBins*sizeof(float));

  

  for (i_bin=min_bin; i_bin<=max_bin; i_bin++){
    PhaseBin[i_bin-min_bin] = (float)i_bin/(float)(NBin-1);
    Profile[i_bin-min_bin] = ProfArray[i_bin];
  }

  cpgbbuf(); 
  cpgeras();
  /* Set world coordinates for current plot */
  cpgswin(x1, x2, FMin(Profile, NewBins, &i_min), 
	  FMax(Profile, NewBins, &i_max));
  cpgline(NewBins, PhaseBin, Profile);
  
  /* Set up plot and labelling preferences */
  cpgbox("BCTSN", 0.0, 0.0, "BCTSVN", 0.0, 0.0);            

  //  if(i_plot == 0){
    sprintf(x_label,"Profile phase");
    sprintf(y_label,"Normalised Flux Density");
    // i_plot++;
    //  }
  cpgmtxt("B", 2.35, 0.5, 0.5, x_label);
  cpgmtxt("L", 3.1, 0.5, 0.5, y_label);
  
  cpgebuf();  
  
  free(PhaseBin);
  free(Profile);


}


/* Construct array out of only good points using wgt array (1s and 0s) */
void GetHistArray(int n_array, float *array, int *wgt, 
		  int *n_hist_array, float *hist_array)
{
  
  int i_array;

  (*n_hist_array)=0;
  
  for (i_array=0; i_array<n_array; i_array++){
    if(wgt[i_array]) hist_array[(*n_hist_array)++] = array[i_array];
    //if((*n_hist_array) == 103) printf("hist_array[102] = %f, array[%d] = %f, wgt = %d, n_array = %d\n",
    //				   hist_array[102], i_array, array[i_array],
    //				      wgt[i_array], n_array);
 }
  

}

  






int GetHist(int n_array, float *array, int n_bin, float *bin_centre, float *hist)
{

  int i_bin, i_array;
  int i_min, i_max, which_bin;
  float bin_size, min_array, max_array;
  
  min_array = FMin(array, n_array, &i_min);
  max_array = FMax(array, n_array, &i_max);
  bin_size = (max_array - min_array)/(float)n_bin;

  //  printf("Made it 1!\n");fflush(stdout);

  /* Initialise histogram by calculating bin centre values and setting counts to zero */
  for (i_bin=0; i_bin<n_bin; i_bin++){
    bin_centre[i_bin] = (double)i_bin*bin_size + 0.5*bin_size + min_array;
    hist[i_bin] = 0.;
  }
  // printf("Made it 2!\n");fflush(stdout);

  // printf("**** n_bin_hist = %d ****\n", n_bin); fflush(stdout);
  /* Drop value into appropriate bin */
  for(i_array=0; i_array<n_array; i_array++){
    //  if(i_array==i_max) // Safeguard so that max value doesn't spill over into next bin
    /* Safeguard so that max value doesn't spill over into next bin */
    /* Worst case, will be = max_array, unless round-off errors cause it to spill 
       slightly over, make condition ">=" */
    if(array[i_array] >= max_array)
      which_bin = n_bin-1;
    else 
      which_bin = (int)(floorf((array[i_array] - min_array) / bin_size));
    // printf("Made it 3!\n");fflush(stdout);
    //printf("array[%d] = %f,  min_array = %f,  bin_size = %f\nWHICH_BIN = %d\n", 
    //	i_array, array[i_array], min_array, bin_size, which_bin);fflush(stdout);
    hist[which_bin] += 1.0;
    //printf("Made it 4!\n");fflush(stdout);
    if(which_bin < 0 || which_bin >= n_bin){
      fprintf(stderr, "Found an out-of-range bin for histogram!\n"); 
      fprintf(stderr, "calculated bin %d, number of bins %d\n", which_bin, n_bin);
      fprintf(stderr, "Value %f, Min = %f, Mx = %f\n", array[i_array], 
	      min_array, max_array);
     return -1;
    }
  }
  
  return 1;

}


void PlotHist(int n_bin_hist, float *bin_centre, float *hist, char *x_label, char *y_label)
{  
  int i_xmin, i_xmax, i_ymax;
  float x1, x2, y1, y2;
  
  x1 = FMin(bin_centre, n_bin_hist, &i_xmin);
  x2 = FMax(bin_centre, n_bin_hist, &i_xmax);
  //  y1 = FMin(hist, n_bin_hist, &i_ymin);
  y1 = 0.;
  y2 = FMax(hist, n_bin_hist, &i_ymax);

  //printf("bin_centre[0] = %f\n", bin_centre[0]);
  // printf("hist[0] = %f\n", hist[0]);
  // printf("n_bin_hist = %d\n", n_bin_hist);
  // printf("x1 = %f, x2 = %f, y1 = %f, y2 = %f\n", x1, x2, y1, y2);

  /* Set world coordinates for current plot */
  cpgswin(x1, x2, y1, y2 + 0.1*(y2-y1));
    

  /* Set up plot and labelling preferences */
  cpgbox("BCTSN", 0.0, 0.0, "BCTSVN", 0.0, 0.0);            
  
  /* Plot line of data array */
  //cpgline(n_bin_hist, bin_centre, hist);
  //MinHist = FMin(ProfRMS, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_min_hist);
  // MaxHist = FMax(ProfRMS, Hdr.redn.RNTimeDumps*Hdr.obs.NChan, &i_max_hist);

  cpgbin(n_bin_hist, bin_centre, hist, 1);
  //  cpghist(Hdr.redn.RNTimeDumps*Hdr.obs.NChan, ProfRMS, MinHist, MaxHist, n_bin_hist, 1);
  
  cpgmtxt("B", 2.35, 0.5, 0.5, x_label);
  cpgmtxt("L", 3.1, 0.5, 0.5, y_label);

}

void PlotAllHist(int n_bin_hist, float *RMSBinVal, float *RMSHist, 
		 float *ShiftBinVal, float *ShiftHist, 
		 float *eShiftBinVal, float *eShiftHist, 
		 float *ScaleBinVal, float *ScaleHist, int Template)
{


  char  x_label[64], y_label[64];
  //  cpgpage();


  /* Open buffer */
  cpgbbuf();  

 
 /* Set up labels for each histogram plot */
  sprintf(x_label,"RMS Flux density");
  sprintf(y_label," ");
  if(Template) cpgpanl(1,1);
  cpgeras();
  PlotHist(n_bin_hist, RMSBinVal, RMSHist, x_label, y_label);

  /* Now do other plots if we are given a template */
  if(Template){
    /* Cross correlation shifts */
    /*    cpgsvp(xaxis_min, xaxis_max, yaxis_max-2.*(yaxis_max-yaxis_min)/4., 
	  -0.03 + yaxis_max-1.*(yaxis_max-yaxis_min)/4.); */
    cpgpanl(1,2);
    sprintf(x_label,"Cross-correlation shift");
    sprintf(y_label," ");
    cpgeras();
    PlotHist(n_bin_hist, ShiftBinVal, ShiftHist, x_label, y_label);
    /* Cross-correlation shift uncertainties */
    /* cpgsvp(xaxis_min, xaxis_max, -0.03 + yaxis_max-3.*(yaxis_max-yaxis_min)/4., 
       -0.06 + yaxis_max-2.*(yaxis_max-yaxis_min)/4.); */
    cpgpanl(1,3);
    sprintf(x_label,"Cross-correlation shift uncertainty");
    sprintf(y_label," ");
    cpgeras();
    PlotHist(n_bin_hist, eShiftBinVal, eShiftHist, x_label, y_label);
    /* Cross-correlation scale factors */
    /*   cpgsvp(xaxis_min, xaxis_max, -0.06 + yaxis_max-4.*(yaxis_max-yaxis_min)/4., 
	 -0.09 + yaxis_max-3.*(yaxis_max-yaxis_min)/4.); */
    cpgpanl(1,4);
    sprintf(x_label,"Cross-correlation scale factor");
    sprintf(y_label," ");
    cpgeras();
    PlotHist(n_bin_hist, ScaleBinVal, ScaleHist, x_label, y_label);

  }
  cpgebuf();

}

/*******************************************************************************
 * Routines to deal with user input for histogram, rms, and profile cuts
 *******************************************************************************/


/* Get and organise input from user for rejection based on histograms */
/* If AcceptInverse is 0, first limit in array can be larger than second, allowing 
   kept region to wrap around right edge of histogram, e.g. for template shift
   values */
void GetHistLimits(float *HistLim, int AcceptWrap)
{
 
  int   i;
  float x_lim_lo, x_lim_hi, temp_limit;
  char  input_str[64];

  fgets(&input_str[0], 64, stdin);
 
  if((strlen(input_str) <= 1 && !strcmp(input_str,"\n"))){
    
    printf("Will leave current limits as before.\n");
  }
  //if(!strncmp(&input_str[0]," ",strlen(input_str))){
  //	printf("All spaces -- Will leave current limits as before.\n");
  //	}
  else{
    for(i=0; i<strlen(input_str); i++){
      if(strncmp(&input_str[i]," ",1)) // i.e. non-space
	break;
    }	
    if(i==strlen(input_str)-1){
      printf("Will leave current limits as before.\n");
      
    }
    else{
      //printf("i = %d\n", i);
      
      
      sscanf(input_str,"%f %f",&x_lim_lo, &x_lim_hi);
      //printf("\n");
      
      if((x_lim_lo==0.0 && x_lim_hi==0.0)){
	printf("Problem with one or more input values. Assuming no change in limits\n");
      }
      else if(x_lim_lo == x_lim_hi){
	printf("Minimum and maximum limits are the same. Will assume there is no ");
	printf("change in limits.\n");
      }
      else{
	/* Quick check that limits are in correct order. If not, switch them, only if 
	   we are not allowing the kept region to wrap around histogram edge */
	if(x_lim_lo > x_lim_hi && !AcceptWrap){
	  temp_limit = x_lim_lo;
	  x_lim_lo = x_lim_hi;
	  x_lim_hi = temp_limit; 
	} 
	/* Now assign limits to array to be returned IF they are inside the current 
	   limits.  Otherwise, just keep them at current limits. */
	if(x_lim_lo > HistLim[0]) 
	  HistLim[0] = x_lim_lo;
	if(x_lim_hi < HistLim[1])
	  HistLim[1] = x_lim_hi;
	
      }
    }
  }
  
}


/* Routine to draw vertical lines to denote where zap limits lie on histogram */
void DrawHistLimits(int n_bin, float *x, float *y, float *lim)
{
  
  int i_dummy; /* Dummy variable to pass into Max/Min function, not needed otherwise */
  float x_min, x_max, y_min, y_max;
    
  x_min = FMin(x, n_bin, &i_dummy);
  x_max = FMax(x, n_bin, &i_dummy);
  y_min = 0.;
  y_max = FMax(y, n_bin, &i_dummy);

  cpgswin(x_min, x_max, 0., y_max);
    //    cpgmove(lim[0], 0.0);
    //    cpgdraw(lim[0], y_max);
  cpgsfs(3);
  cpgshs(45.0, 2.0, 0.0);
  if(lim[0] > lim[1]){  /* region being kept is wrapping around edge of histogram */
    cpgrect(lim[1], lim[0], y_min-0.1, y_max+0.1);
  }
  else{
    cpgrect(x_min-0.1, lim[0], y_min-0.1, y_max+0.1);
    
    //    cpgmove(lim[1], 0.0);
    //    cpgdraw(lim[1], y_max);
    cpgshs(-45.0, 2.0, 0.0);
    cpgrect(lim[1], x_max+0.1, y_min-0.1, y_max+0.1);
  }
}

/* Simple routine to get a yes or no input from user */
int GetYesNo(char *Question)
{
 
  char input_str[64];
  int retval=-1;

  while(retval != YES && retval != NO){
    printf("%s (y/n) [y]: ", Question);
    fgets(&input_str[0], 64, stdin);
    
    if(!strncmp(input_str,"\n", 1) || 
       //       !strncmp(input_str, "y\n", 2) || 
       !strncmp(input_str, "y", 1) || 
       //       !strncmp(input_str, "yes\n", 4)){
       !strncmp(input_str, "yes", 3)){
      //printf("You said YES\n");
      retval = YES;
    }
    //    else if(!strncmp(input_str, "n\n", 2) || 
    //	    !strncmp(input_str, "no\n", 3)){
    else if(!strncmp(input_str, "n", 1) || 
	    !strncmp(input_str, "no", 2)){
      //printf("You said NO\n");
      retval = NO;
    }
    else{
      printf("Unrecognized response. Try again... \n");
    }
  }

  return retval;

}


void PlotAllProf(float *ProfChan, int *ProfChanMask, 
		 float *ProfDump, int *ProfDumpMask,
		 struct ASPHdr Hdr, 
		 float xaxis_min, float xaxis_max, float yaxis_min, float yaxis_max,
		 struct gray *ProfChanGray, struct plot_info *ProfChanInfo,
		 struct gray *ProfDumpGray, struct plot_info *ProfDumpInfo)
{

  float x1, x2, y1, y2;

  cpgeras();

  /* Left-hand panel */
  x1 = xaxis_min;
  x2 = (xaxis_min + xaxis_max)/2. - 0.06;
  y1 = yaxis_min;
  y2 = yaxis_max;

  /*  printf("**** Checking ProfChan *****\n");
  for (i_chan=0; i_chan<Hdr.obs.NChan; i_chan++) {
    for (i=0; i<Hdr.redn.RNBinTimeDump; i++) {
    }
    } 
    printf("\n"); */

  /* Add in a half-channel at each end to ensure that labels refer to 
     channel centres */
  GetGray(ProfChan, ProfChanMask, Hdr.redn.RNBinTimeDump*Hdr.obs.NChan, 
	  1, Hdr.redn.RNBinTimeDump, 1, Hdr.obs.NChan, 
	  0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
	  Hdr.obs.ChanFreq[0], // - Sideband*0.5*ChanBW, 
	  Hdr.obs.ChanFreq[Hdr.obs.NChan-1], // + Sideband*0.5*ChanBW,
	  ProfChanGray, ProfChanInfo);

  sprintf(ProfChanInfo->x_label,"Profile phase");
  sprintf(ProfChanInfo->y_label,"Observing frequency (MHz)");
  
  //printf("z = [%f, %f]", ProfChanGray->z[0], ProfChanGray->z[1]);
  cpgsvp(x1, x2, y1, y2); // half the window

  //cpgeras();

  PlotGray(ProfChanGray, ProfChan, ProfChanInfo);

  /* Now do the other profile grayscale, on the right side */

  /* Right-hand panel */
  x1 = (xaxis_min + xaxis_max)/2. + 0.06;
  x2 = xaxis_max;
  y1 = yaxis_min;
  y2 = yaxis_max;

  GetGray(ProfDump, ProfDumpMask, Hdr.redn.RNBinTimeDump*Hdr.redn.RNTimeDumps, 
	  1, Hdr.redn.RNBinTimeDump, 1, Hdr.redn.RNTimeDumps, 
	  0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
	  0.0, (float)(Hdr.redn.RNTimeDumps-1),
	  ProfDumpGray, ProfDumpInfo);
  
  sprintf(ProfDumpInfo->x_label,"Profile phase");
  sprintf(ProfDumpInfo->y_label,"Integration number");

  cpgsvp(x1, x2, y1, y2); // half the window
  
  PlotGray(ProfDumpGray, ProfDump, ProfDumpInfo);

  return;

}


int ProfZoom(float *ProfArray, int *ProfMask, struct ASPHdr Hdr, 
	     float *ZoomVal, int *ZoomIndex, int prof_mode)
{

  int   i_array, StartInd, EndInd;
  float ZoomDiff[2], temp_diff, StartVal, EndVal;
  struct gray ZoomGray;
  struct plot_info ZoomInfo;

  /* At this point, ZoomVal[1] > ZoomVal[0] */

  /* Find out what the zoom limits are and find out to which array elements 
     they correspond */
  
  
  /* If we are doing the frequency profiles, we need to make sure of which channel the click 
     is within */
  if (prof_mode==FREQMODE){
    //ZoomIndex[0] = -9999;
    //ZoomIndex[1] = -9999;
    ZoomDiff[0] = fabsf(ZoomVal[0] - Hdr.obs.ChanFreq[0]);
    ZoomDiff[1] = fabsf(ZoomVal[1] - Hdr.obs.ChanFreq[0]);
    /* Start at index 1 since we've already calculated difference for the 0th element */
    for(i_array=1; i_array<Hdr.obs.NChan; i_array++){
      if( (temp_diff = fabsf(ZoomVal[0] - Hdr.obs.ChanFreq[i_array])) < ZoomDiff[0]){ 
	ZoomIndex[0] = i_array;
	ZoomDiff[0] = temp_diff;
      }
      if( (temp_diff = fabsf(ZoomVal[1] - Hdr.obs.ChanFreq[i_array])) < ZoomDiff[1]){
	ZoomIndex[1] = i_array;     
	ZoomDiff[1] = temp_diff;
      }
    }
    if(ZoomIndex[0] < 0 || ZoomIndex[1] < 0){
      fprintf(stderr, "Error: Could not calculate clicked channels.\n");
      return -1;
    }
    else{
      ZoomVal[0] = Hdr.obs.ChanFreq[ZoomIndex[0]];
      ZoomVal[1] = Hdr.obs.ChanFreq[ZoomIndex[1]];

      printf("You have chosen to zoom between channels %d (%.3f) and %d (%.3f).\n",
	     ZoomIndex[0], Hdr.obs.ChanFreq[ZoomIndex[0]], 
	     ZoomIndex[1], Hdr.obs.ChanFreq[ZoomIndex[1]]);

      /* Now plot zoomed grayscale plot */
      sprintf(ZoomInfo.x_label,"Profile phase");
      sprintf(ZoomInfo.y_label,"Observing frequency (MHz)");

    }
  }

  /* If we are in time-added profile mode */
  else if(prof_mode==TIMEMODE){

    /* Test that Chosen Zoom value clicks are within the current plotted values.  
       If not, keep extreme values as they are */
    if(ZoomVal[0] > (float)ZoomIndex[0])
      ZoomIndex[0] = (int)(ZoomVal[0] + 0.5);
    
    if(ZoomVal[1] < (float)ZoomIndex[1])
      ZoomIndex[1] = (int)(ZoomVal[1] + 0.5);
    
    ZoomVal[0] = (float)ZoomIndex[0];
    ZoomVal[1] = (float)ZoomIndex[1];

    printf("You have chosen to zoom between subints %d and %d.\n",
	   ZoomIndex[0],  ZoomIndex[1]);
   
    /* Now plot zoomed grayscale plot */
    sprintf(ZoomInfo.x_label,"Profile phase");
    sprintf(ZoomInfo.y_label,"Integration number");
    
    
  }
  else{
    fprintf(stderr, "Error: Unrecognized profile mode.\n");
    return -1;
  }

  cpgeras();
  /* Since it is a simple range, we can easily send over starting pointer and array length 
     to the GetGray routine, and get back our grayscale arrays (which are pre-allocated) */
  if(ZoomIndex[0] < ZoomIndex[1]){
    StartInd = ZoomIndex[0];
    EndInd = ZoomIndex[1];
    StartVal = ZoomVal[0];
    EndVal = ZoomVal[1];
  }
  else{
    StartInd = ZoomIndex[1];
    EndInd = ZoomIndex[0];	
    StartVal = ZoomVal[1];
    EndVal = ZoomVal[0];
  }
  
  
  
  GetGray(&ProfArray[StartInd*Hdr.redn.RNBinTimeDump], 
	  &ProfMask[StartInd*Hdr.redn.RNBinTimeDump], 
	  (EndInd-StartInd+1)*Hdr.redn.RNBinTimeDump, 
	  1, Hdr.redn.RNBinTimeDump, 1, 1+(EndInd-StartInd), 
	  0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
	  StartVal, EndVal,
	  &ZoomGray, &ZoomInfo);
  PlotGray(&ZoomGray, &ProfArray[StartInd*Hdr.redn.RNBinTimeDump], &ZoomInfo);
  
  return 1;
  
}


int ProfZap(float *ProfArray, int *ProfMask, struct ASPHdr Hdr, 
	    float *ZapVal, int *ZoomIndex, int *ZapIndex, int *ZapFlag, int prof_mode)
{

  int    i_array, i_mask, i_zap, StartInd, EndInd;
  int    *TempMask, temp_ind;
  float  ZoomVal[2], ZapDiff[2], temp_diff, StartVal, EndVal;
  struct gray Gray;
  struct plot_info Info;
  

  if (prof_mode==FREQMODE){
    ZoomVal[0] = Hdr.obs.ChanFreq[ZoomIndex[0]];
    ZoomVal[1] = Hdr.obs.ChanFreq[ZoomIndex[1]];
    
    TempMask = (int *)malloc(Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));
    memcpy(TempMask, ProfMask, Hdr.obs.NChan*Hdr.redn.RNBinTimeDump*sizeof(int));


    ZapDiff[0] = fabsf(ZapVal[0] - Hdr.obs.ChanFreq[0]);
    ZapDiff[1] = fabsf(ZapVal[1] - Hdr.obs.ChanFreq[0]);
    /* Start at index 1 since we've already calculated difference for the 0th element */
    for(i_array=1; i_array<Hdr.obs.NChan; i_array++){
      if( (temp_diff = fabsf(ZapVal[0] - Hdr.obs.ChanFreq[i_array])) < ZapDiff[0]){ 
	ZapIndex[0] = i_array;
	ZapDiff[0] = temp_diff;
      }
      if( (temp_diff = fabsf(ZapVal[1] - Hdr.obs.ChanFreq[i_array])) < ZapDiff[1]){
	ZapIndex[1] = i_array;     
	ZapDiff[1] = temp_diff;
      }
    }
    if(ZapIndex[0] < 0 || ZapIndex[1] < 0){
      fprintf(stderr, "Error: Could not calculate clicked channels.\n");
      return -1;
    }
    else{
      ZapVal[0] = Hdr.obs.ChanFreq[ZapIndex[0]];
      ZapVal[1] = Hdr.obs.ChanFreq[ZapIndex[1]];

      printf("You have chosen to zap channels within the range %d (%.3f) and %d (%.3f).\n",
	     ZapIndex[0], Hdr.obs.ChanFreq[ZapIndex[0]], 
	     ZapIndex[1], Hdr.obs.ChanFreq[ZapIndex[1]]);

      /* Now plot zoomed grayscale plot */
      sprintf(Info.x_label,"Profile phase");
      sprintf(Info.y_label,"Observing frequency (MHz)");

    }
    
  }
  /* If we are in time-added profile mode */
  else if(prof_mode==TIMEMODE){
    ZoomVal[0] = (float)ZoomIndex[0];
    ZoomVal[1] = (float)ZoomIndex[1];
   
    TempMask = (int *)malloc(Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));
    memcpy(TempMask, ProfMask, Hdr.redn.RNTimeDumps*Hdr.redn.RNBinTimeDump*sizeof(int));

    /* Test that Chosen Zoom value clicks are within the current plotted values.  
       If not, keep extreme values as they are */
    if(ZoomIndex[0] < ZoomIndex[1]) {
      for (i_zap=0; i_zap<2; i_zap++){
	if(ZapVal[i_zap] > (float)ZoomIndex[0] && ZapVal[i_zap] < (float)ZoomIndex[1]) {
	  ZapIndex[i_zap] = (int)(ZapVal[i_zap] + 0.5);
	}
	else{
	  if(ZapVal[i_zap] < (float)ZoomIndex[0]) ZapIndex[i_zap] = ZoomIndex[0];
	  else ZapIndex[i_zap] = ZoomIndex[1];
	}
      }
    }
    else{
      for (i_zap=0; i_zap<2; i_zap++){
	if(ZapVal[i_zap] < (float)ZoomIndex[0] && ZapVal[i_zap] > (float)ZoomIndex[1]) {
	  ZapIndex[i_zap] = (int)(ZapVal[i_zap] + 0.5);
	}
	else{
	  if(ZapVal[i_zap] > (float)ZoomIndex[0]) ZapIndex[i_zap] = ZoomIndex[0];
	  else ZapIndex[i_zap] = ZoomIndex[1];
	}
      }
      
    }

   
    ZapVal[0] = (float)ZapIndex[0];
    ZapVal[1] = (float)ZapIndex[1];

    printf("You have chosen to zoom between subints %d and %d.\n",
	   ZapIndex[0],  ZapIndex[1]);
   
    /* Now plot zoomed grayscale plot */
    sprintf(Info.x_label,"Profile phase");
    sprintf(Info.y_label,"Integration number");
    
  
  }
  else{
    fprintf(stderr, "Error: Unrecognized profile mode.\n");
    return -1;
  }
 
  /* Now get Zoom and Zap Indices int he right order */
  if(ZoomIndex[0] < ZoomIndex[1]){
    StartInd = ZoomIndex[0];
    EndInd = ZoomIndex[1];
    StartVal = ZoomVal[0];
    EndVal = ZoomVal[1];
  }
  else{
    StartInd = ZoomIndex[1];
    EndInd = ZoomIndex[0];	
    StartVal = ZoomVal[1];
    EndVal = ZoomVal[0];
  }

  if(ZapIndex[0] > ZapIndex[1]){
    temp_ind = ZapIndex[0];
    ZapIndex[0] = ZapIndex[1];
    ZapIndex[1] = temp_ind;
  }
  /* Now set ZapFlags for the range of frequencies/subints here */
  for (i_zap=ZapIndex[0]; i_zap<=ZapIndex[1]; i_zap++) 
    ZapFlag[i_zap]=1;

  /* Now add these to a temporary Mask Array based on the current Mask, so it is not 
     yet made permanent */
  for (i_mask=ZapIndex[0]*Hdr.redn.RNBinTimeDump; 
       i_mask<(ZapIndex[1]+1)*Hdr.redn.RNBinTimeDump; i_mask++){
    ProfMask[i_mask] = 0;
  }
  
  /* Now plot temporarily-masked zapped plot */
  cpgeras();
  GetGray(&ProfArray[StartInd*Hdr.redn.RNBinTimeDump], 
	  &ProfMask[StartInd*Hdr.redn.RNBinTimeDump], 
	  (EndInd-StartInd+1)*Hdr.redn.RNBinTimeDump, 
	  1, Hdr.redn.RNBinTimeDump, 1, 1+(EndInd-StartInd), 
	  0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
	  StartVal, EndVal,
	  &Gray, &Info);
  PlotGray(&Gray, &ProfArray[StartInd*Hdr.redn.RNBinTimeDump], &Info);
  
  printf("ZAPINDEX = %d  %d\n\n", ZapIndex[0], ZapIndex[1]);

  free(TempMask);
  return 1;
}

int UndoZap(float *ProfArray, int *ProfMask, struct ASPHdr Hdr, 
	    int *ZoomIndex, int *ZapIndex, int *ZapFlag, int prof_mode)
{

  int    i_mask, i_zap, StartInd, EndInd;
  float  ZoomVal[2], StartVal, EndVal;
  struct gray Gray;
  struct plot_info Info;

  /* Change Mask values between Zap Indices back to 1 */
  for (i_mask=ZapIndex[0]*Hdr.redn.RNBinTimeDump; 
       i_mask<(ZapIndex[1]+1)*Hdr.redn.RNBinTimeDump; i_mask++){
    ProfMask[i_mask] = 1;
  }
  /* Change flag values between Zap Indices back to 0 */
  for (i_zap=ZapIndex[0]; i_zap<=ZapIndex[1]; i_zap++)
    ZapFlag[i_zap]=0;

  
  /* Now re-plot grayscale plot in current zoom */

  if (prof_mode==FREQMODE){
    ZoomVal[0] = Hdr.obs.ChanFreq[ZoomIndex[0]];
    ZoomVal[1] = Hdr.obs.ChanFreq[ZoomIndex[1]];

    sprintf(Info.x_label,"Profile phase");
    sprintf(Info.y_label,"Observing frequency (MHz)");
  }
  /* If we are in time-added profile mode */
  else if(prof_mode==TIMEMODE){
    ZoomVal[0] = (float)ZoomIndex[0];
    ZoomVal[1] = (float)ZoomIndex[1];

    sprintf(Info.x_label,"Profile phase");
    sprintf(Info.y_label,"Integration number");
  }
  else{
    fprintf(stderr, "Error: Unrecognized profile mode.\n");
    return -1;
  }


  /* Now get Zoom and Zap Indices in the right order */
  if(ZoomIndex[0] < ZoomIndex[1]){
    StartInd = ZoomIndex[0];
    EndInd = ZoomIndex[1];
    StartVal = ZoomVal[0];
    EndVal = ZoomVal[1];
  }
  else{
    StartInd = ZoomIndex[1];
    EndInd = ZoomIndex[0];	
    StartVal = ZoomVal[1];
    EndVal = ZoomVal[0];
  }

  /* Now plot */
  cpgeras();
  GetGray(&ProfArray[StartInd*Hdr.redn.RNBinTimeDump], 
	  &ProfMask[StartInd*Hdr.redn.RNBinTimeDump], 
	  (EndInd-StartInd+1)*Hdr.redn.RNBinTimeDump, 
	  1, Hdr.redn.RNBinTimeDump, 1, 1+(EndInd-StartInd), 
	  0.,  (float)(Hdr.redn.RNBinTimeDump-1), 
	  StartVal, EndVal,
	  &Gray, &Info);
  PlotGray(&Gray, &ProfArray[StartInd*Hdr.redn.RNBinTimeDump], &Info);
 
  printf("UNZAPINDEX = %d  %d\n\n", ZapIndex[0], ZapIndex[1]);

  return 1;

}

int ReadZap(struct ASPHdr *hdr, char *ZapIn, int *ProfWgt){

  int    i_chan, i_dump, i_array;
  double freq; 
  char   ZapLine[128];
  FILE   *FZapIn;

  /* Open the input zap file */
  if((FZapIn = fopen(ZapIn, "r")) == NULL){
    fprintf(stderr, "Error in opening output zap file %s.\n", ZapIn);
    exit(1);
  }

  while (fgets(ZapLine, 128, FZapIn) != NULL){
    sscanf(ZapLine, "%d %lf", &i_dump, &freq);
    /* Check first whether frequency is negative.  If so, all channel of the 
       given dump number are flagged for zapping */
    if(freq < 0.0){
      for(i_chan=0; i_chan<hdr->obs.NChan; i_chan++){
	i_array = i_chan*hdr->redn.RNTimeDumps + i_dump;
	ProfWgt[i_array] = 0;
      }
    }
    else {
      /* Now we know frequency is not negative. So, map frequency to integer 
	 channel array number */
      if((i_chan = Freq2Chan(freq, hdr->obs.ChanFreq, hdr->obs.NChan)) < 0) {
	/* If mapping unsuccessful, quit out with error */
	fprintf(stderr, "Error in scan omission file. Frequency %f does not ",
		freq);
	fprintf(stderr, "map to any channel for this data file. \n");
	return -1;
      }
      
      /* Now check whether dump number is negative. If so, all dumps at the 
	 given frequency are flagged for omission */
      if(i_dump < 0){
	for (i_dump=0; i_dump<hdr->redn.RNTimeDumps; i_dump++){
	  i_array = i_chan*hdr->redn.RNTimeDumps + i_dump;
	  ProfWgt[i_array] = 0;
	}	
      }
      else{
	/* Finally, if we got to here, then neither frequency nor dump number 
	   are negative.  In this case, flag this dump/channel combination 
	   as one for omission */
	i_array = i_chan*hdr->redn.RNTimeDumps + i_dump;
	ProfWgt[i_array] = 0;
      }
      
    }
  }
  
  fclose(FZapIn);
  
  return 1;

}
