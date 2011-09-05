/* ========================================================================= */
/* GetPhases                                                                 */
/*     Will take in a (for now) ascii profile in asp format and some input   */
/*     options, and will output a pulse width and pulse width error.         */
/*                                                                           */
/*     User has option to plot the input profile and show the width lines    */
/*     overplotted.                                                          */
/*                                                                           */
/* R. Ferdman, 2011 June 14                                                  */
/* ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "ASPDefs.h"
#include "WidthCmdLine.h"
#include "ASPCommon.h"
#include "cpgplot.h"
#include "cpgplot_def.h"

int NormProf(int, float *, char *);
// double XInterp(int, float *, int, int, double, double, double *);
double XInterp2D(double *, double *, double, double, double *);
void PlotProf(int, float *, float, float);

int main(int argc, char **argv)
{

  int i_file, i_prof, i_bin, i_pk, i_range;
  int i_low, i_high, i_ymin, i_ymax, n_range;
  int bin[NBINMAX], bin_seg[2], NBins, TempNBins;
  int at_height;
  double iMJDStart, sMJDStart;
  double Duty, PeakHeight, Height, SBase;
  double FinalMask[NBINMAX];
  double ProfLow, ProfHigh;
  double Width, WidthErr;
  double phase_seg[2];
  /* bounding valus going into estimating interpolated x value: */
  double x_int[2], y_int[2]; 
  double x0, x1, x0_err, x1_err;
  double y_min, y_max;

  int dev_prof_plot, l_width;
  float xaxis_min, xaxis_max, yaxis_min, yaxis_max;

  char Headerline[256], PSRName[12], outfile[256];
  struct StdProfs Profile;
  Cmdline *Cmd;

  FILE **fpout;

  int ijunk;
  double fjunk;
  char sjunk[16];

  /* Get command line variables */
  Cmd = parseCmdline(argc, argv);  

  /* Normally use this somewhere, and not showOptionValues */
  Cmd->tool = Cmd->tool;

  /* Begin plottingsetting up pgplot section */

  if ((dev_prof_plot = cpgopen("/xs")) < 1) {
    printf("Warning:  pgplot device could not be initialized.\n");
    exit(1);
  }
  cpgpap(8., 0.618034);
  cpgsch(1.8);  

  /* Set limiting x and y-axis coordinates for all the windows */
  xaxis_min = 0.11;
  xaxis_max = 0.97;
  yaxis_min = 0.12;
  yaxis_max = 0.97;

  /* End pgplot setup */

  /* Determine the number of phase ranges we will have */
  n_range = Cmd->PhaseSplitC + 1;
  fpout = (FILE **)malloc(n_range*sizeof(FILE *));
  /* Open output file(s) for each range */
  if(n_range == 1){
    sprintf(outfile, "%s.w%.0lf.dat", Cmd->PSRName, Cmd->PercentHeight);
    if ((fpout[0] = fopen(outfile,"w"))==NULL){
      fprintf(stderr, "Could not open file %s.  Exiting.\n", outfile);
      exit(2);
    }
  }
  else{
    for (i_range=0; i_range<n_range; i_range++){
      if(i_range==0) 
	phase_seg[0] = 0.0;
      else 
	phase_seg[0] = Cmd->PhaseSplit[i_range-1];
      if(i_range==n_range-1) 
	phase_seg[1] = 1.0;
      else 
	phase_seg[1] = Cmd->PhaseSplit[i_range];
      sprintf(outfile, "%s.w%.0lf.%.1f_%.1f.dat", Cmd->PSRName, 
	      Cmd->PercentHeight, phase_seg[0], phase_seg[1]);
      if ((fpout[i_range] = fopen(outfile,"w"))==NULL){
	fprintf(stderr, "Could not open file %s.  Exiting.\n", outfile);
	exit(2);
      }
    }
  }

  /* Will have to change this when we use fits files, with multiple profiles 
     per file (and thus can't determine size of Profile by number of input 
     data files */

  for (i_prof=0; i_prof<Cmd->InfileC; i_prof++){

    if(Cmd->VerboseP)
      printf("Opening file %s\n", Cmd->Infile[i_prof]);
    

    if ( ReadASPAsc(Cmd->Infile[i_prof], &Headerline[0], bin,  
		    &Profile, &NBins) < 0) {
      printf("Error in reading file %s.\n",Cmd->Infile[i_prof]);
      fflush(stdout);
      exit(1);
    }

    /* Read in header variables of use to this routine:  MJD and PSR name */
    sscanf(Headerline, "%s%lf%lf%lf%d%lf%lf%d%d%d%s%lf",
	   sjunk, &iMJDStart, &sMJDStart, &fjunk, &ijunk, &fjunk, &fjunk, 
	   &ijunk, &ijunk, &ijunk, PSRName, &fjunk);

    /* Normalise profile */
    if(NormProf(NBins, Profile.rstds, Cmd->PSRName) < 0){
      printf("Profile could not be normalised.  Have left as is.\n");
    }
    
    /* Plot profile */
    cpgslct(dev_prof_plot);
    cpgpage();
    cpgsvp(xaxis_min, xaxis_max, yaxis_min, yaxis_max);
    
    PlotProf(NBins, Profile.rstds, 0.0, 1.0);
    /* End plot profile */

    /* Get off-pulse rms and also get off-pulse phase bins using BMask */
    Duty = DutyLookup(Cmd->PSRName);
    BMask(Profile.rstds,&NBins,&Duty,FinalMask);
    Baseline(Profile.rstds,FinalMask,&NBins,&SBase,&Profile.Srms);
    Profile.Srms = 1.0/Profile.Srms;


    /* Now split off into separate phase ranges, and get widths for component 
       of profile within each of these ranges.  If n_range == 1, then use 
       entire range. */
    
    for (i_range=0; i_range<n_range; i_range++){
      
      /* Assign phase bin edges (in bins, not phase), keeping in mind the 
	 0 and 1 edges for the first and last bin, respectively */
      if(i_range==0) 
	bin_seg[0] = 0;
      else
	bin_seg[0] = (int)((double)(NBins-1)*Cmd->PhaseSplit[i_range-1]);
      if(i_range==n_range-1)
	bin_seg[1] = NBins-1;
      else
	bin_seg[1] = (int)((double)(NBins-1)*Cmd->PhaseSplit[i_range]);
      
      TempNBins = bin_seg[1] - bin_seg[0] + 1;
      
      //      PeakHeight =  FindPeak(Profile.rstds,&NBins,&i_pk);
      PeakHeight =  FindPeak(&Profile.rstds[bin_seg[0]],&TempNBins,&i_pk);
      Profile.SNR = PeakHeight*Profile.Srms;
      Height = 0.01*Cmd->PercentHeight*PeakHeight;
    
      if (Cmd->VerboseP){
	printf("Requested width: %lf\nPeak Height: %lf\nWidth Height: %lf\nRMS: %lf\n\n",
	       Cmd->PercentHeight, PeakHeight, Height, Profile.Srms);
      }


      /* Need to find phases at which profile is closest to requested
	 pulse height */
    
      /* Do this by running two points, one from each side, along the 
	 profile until one point is just below and one is just above the 
	 requested fractional pulse height, then interpolate between the two */
       
      i_bin=bin_seg[0];
      at_height=0;
      /* Now start from the left side, going right */
      while(!at_height){
	ProfLow = Profile.rstds[i_bin] - Height;
	ProfHigh = Profile.rstds[i_bin+1] - Height;
	if(ProfLow<=0. && ProfHigh>=0.){
	  at_height=1;
	  i_low = i_bin;
	  i_high = i_bin+1;
	}
	else{
	  i_bin++;
	}
	//	if(i_bin == NBins) {
	if(i_bin == bin_seg[1]) {
	  fprintf(stderr, "Could not get left side of profile width.\n");
	  exit(2);
	}
      }
      /* Now interpolate between the two bins to get the phase at the 
	 requested pulse height */
      //      x0 = XInterp(NBins, Profile.rstds, i_low, i_high, Height, Profile.Srms, &x0_err);
      x_int[0] = (double)i_low;
      x_int[1] = (double)i_high;
      y_int[0] = (double)Profile.rstds[i_low];
      y_int[1] = (double)Profile.rstds[i_high];
      x0 = XInterp2D(x_int, y_int, Height, Profile.Srms, &x0_err);
      printf("Pre-interp phases:  %lf  ", (double)i_bin/(double)NBins);

      /* Now start from the right side, going left */
      //      i_bin=NBins-1;
      i_bin=bin_seg[1];
      at_height=0;
      while(!at_height){
	ProfLow = Profile.rstds[i_bin] - Height;
	ProfHigh = Profile.rstds[i_bin-1] - Height;
	if(ProfLow<=0. && ProfHigh>=0.){
	  at_height=1;
	  i_low = i_bin;
	  i_high = i_bin-1;
	}
	else{
	  i_bin--;
	}
	if(i_bin < bin_seg[0]) {
	  fprintf(stderr, "Could not get right side of profile width.\n");
	  exit(2);
	}
      }
      /* Now interpolate between the two bins to get the phase at the 
	 requested pulse height */
      //      x1 = XInterp(NBins, Profile.rstds, i_low, i_high, Height, Profile.Srms, &x1_err);
      x_int[0] = (double)i_low;
      x_int[1] = (double)i_high;
      y_int[0] = (double)Profile.rstds[i_low];
      y_int[1] = (double)Profile.rstds[i_high];
      x1 = XInterp2D(x_int, y_int, Height, Profile.Srms, &x1_err);
      printf("%lf\n\n", (double)i_bin/(double)NBins);
    

      /* Width is just the difference between the two */
      if(x1 > x0) {
	Width = (x1 - x0)/((double)NBins);
      }
      else {   /* Need to deal with case of profile wrapping through phase 1 and 
		  back from 0 again */
	Width = (x1 - x0 + (double)NBins)/((double)NBins);
      }
	WidthErr = sqrt(x0_err*x0_err + x1_err*x1_err)/((double)NBins);


      printf("%s   %.15lf  %s   %lf   %lf\n", Cmd->Infile[i_prof], 
	     iMJDStart+(sMJDStart/86400.0), Cmd->PSRName, Width, WidthErr);
      fprintf(fpout[i_range], "%.15lf  %lf   %lf\n", 
	      iMJDStart+(sMJDStart/86400.0), Width, WidthErr);

      /* Draw vertical lines at phases between which width will be calculated */
      cpgsci(c_red);
      cpgsls(l_dashed);
      cpgqlw(&l_width); /* Get current line width */
      cpgslw(5*l_width);
    
      y_min = FMin(Profile.rstds, NBins, &i_ymin);
      y_max = FMax(Profile.rstds, NBins, &i_ymax);

      cpgmove((float)x0/(float)NBins, y_min-0.1);
      cpgdraw(x0/(float)NBins, Height);
      cpgmove((float)x1/(float)NBins, y_min-0.1);
      cpgdraw(x1/(float)NBins, Height);
      if(x1 > x0){
	cpgmove((float)x0/(float)NBins, Height);
	cpgdraw(x1/(float)NBins, Height);
      }
      else {
	cpgmove((float)x0/(float)NBins, Height);
	cpgdraw(1.01, Height);
	cpgmove((float)x1/(float)NBins, Height);
	cpgdraw(-0.01, Height);     
      }
      /* Faintly draw dotted line where phase splitting occurs */
      /* Comapring i_range to n_range and not n_range-1 here since 
	 i_range will have been already augmented by 1 here */      
      if(i_range<n_range-1){ 
	cpgsci(c_blue);
	cpgsls(l_dotted);
	cpgslw(3*l_width);
	cpgmove((float)Cmd->PhaseSplit[i_range], y_min-0.1);
	cpgdraw((float)Cmd->PhaseSplit[i_range], y_max+0.1);
      }

      cpgsci(1); /* Back to default */
      cpgsls(l_full);
      cpgslw(l_width); /* reset to default line width) */
    }
    

    /* End plotting section */
  }
  
  for (i_range=0; i_range<n_range; i_range++)
    fclose(fpout[i_range]);
 
  exit(0);

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

#if 0
double XInterp(int NBins, float *Profile, int i_low, int i_high, double y_val, 
	       double rms, double *x_err)
{

  double x_interp, slope;

  slope = (Profile[i_high] - Profile[i_low])/((double)(i_high - i_low));

  x_interp = (double)i_low + (y_val - Profile[i_low])/slope;

  *x_err = rms*sqrt( (1.0/(slope*slope)) 
		    + 2.0*(y_val-Profile[i_low])*(y_val-Profile[i_low])
		    / (slope*slope*slope*slope) );
    
  return x_interp;
}
#endif


double XInterp2D(double *x, double *y, double y_val, 
	       double y_err, double *x_err)
{

  double x_interp, slope;

  slope = (y[1] - y[0])/(x[1] - x[0]);

  x_interp = x[0] + (y_val - y[0])/slope;

  *x_err = y_err*sqrt( (1.0/(slope*slope)) 
		    + 2.0*(y_val-y[0])*(y_val-y[0])
		    / (slope*slope*slope*slope) );
    
  return x_interp;

}
