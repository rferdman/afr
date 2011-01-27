/* .h files, define's, and structures used by all programs */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "ASPDefs.h"
#include "ASPHeader.h"
#include "fitsio.h"
#include "polyco.h"



struct RunVars {
  int    NBins;
  int    NBinRead;
  int    NDumps;
  int    BinDown;
  int    NBinsOut;
  int    MakeRaw;
  int    MakeStd;
  int    FlipPA;
  int    NoBase;
  int    Scale;
  int    Header;
  int    Verbose;
  int    Swap;
  int    Dedisp;
  int    OldFits;
  int    Cal;
  int    ThetaBBFlag;
  int    NScanOmit;
  int    NumAllDumpOmit;
  int    AddChans;
  int    AddDumps;  
  int    SetsofChans2Add;
  int    Sideband;
  int    NumEffChans;  
  int    NumEffDumps;
  int    NOutDumps;
  int    NOutChans;
  int    *FirstChanAdd;
  //  int    FirstChanAdd[NCHMAX];
  int    *LastChanAdd;
  //  int    LastChanAdd[NCHMAX];
  int    *MinChans2Add;
  //  int    MinChans2Add[NCHMAX];
  int    *MaxChans2Add;
  //  int    MaxChans2Add[NCHMAX];
  int    *NumChans2Add;
  //  int    NumChans2Add[NCHMAX];
  int    *CurZapChan;
  //  int    CurZapChan[NCHMAX];
  int    *ZapChan;
  //  int    ZapChan[NCHMAX];
  int    *DumpOmit;
  //  int    DumpOmit[MAXOMIT];
  int    *ChanOmit;
  //  int    ChanOmit[MAXOMIT];
  int    *TotScans;
  //int    TotScans[MAXDUMPS*NCHMAX];
  int    *TotOmit;
  //int    TotOmit[MAXDUMPS*NCHMAX];
  int    *AllDumpOmit;
  //  int    AllDumpOmit[MAXOMIT];
  int    *OmitFlag;
  //int    OmitFlag[MAXDUMPS*NCHMAX];
  float  *MM;
  //  float  MM[16*NCHMAX];
  double FSky;
  double *ThetaBB;
  //  double ThetaBB[NCHMAX];
/*   double CalValue[2]; */
  char   Infile[256];
  char   OutfileRoot[256];
  char   Stdfile[256];
  char   ThetaBBfile[256];
  char   Muellerfile[256];
  char   Source[12];
};

struct CalVars {
  int    CalStrength;;
  int    CalMethod;
  int    NOmit;
  int    NumCalChans;
  int    CalIndex[NCHMAX];
  double Flux;
  double Tcal[2];
  double JyPerCal[2];
  double Gain;
  double GainOnPulsar;
  double AddAngle;
  double ForcePhase;
  double FOmit[NCHMAX];
  char   Calfile[256];
  char   ContCalfile[256];
};

struct StdProfs {
  float  rstds[NBINMAX];
  float  rstdq[NBINMAX];
  float  rstdu[NBINMAX];
  float  rstdv[NBINMAX];
  //  float  rstdinv[NBINMAX];
  float  stds[NBINMAX][2];
  float  stdq[NBINMAX][2];
  float  stdu[NBINMAX][2];
  float  stdv[NBINMAX][2];
  //  float  stdinv[NBINMAX][2];
  float  stdlin[NBINMAX];
  float  stdphi[NBINMAX];
  float  stdphierr[NBINMAX];
  //  float  stdamp[NBINMAX];
  //  float  stdpha[NBINMAX];
  //  float  stdampinv[NBINMAX];
  //  float  stdphainv[NBINMAX];
  float  stdamps[NBINMAX];
  float  stdampq[NBINMAX];
  float  stdampu[NBINMAX];
  float  stdampv[NBINMAX];
  float  stdphas[NBINMAX];
  float  stdphaq[NBINMAX];
  float  stdphau[NBINMAX];
  float  stdphav[NBINMAX];
  //  float  profmax[NBINMAX];
  //  float  profmin[NBINMAX];
  //  float  stdampx[NBINMAX];
  //  float  stdampn[NBINMAX];
  //  float  stdphax[NBINMAX];
  //  float  stdphan[NBINMAX];
  double SNR;
  double Srms;
  double Qrms;
  double Urms;
  double Vrms;
  double Linrms;
};

struct Telescope {
  double Lat;
  double Long;
};

/* Functions */
void   cprofc(float *, int, float *, float *);
void   uncprofc(float *, float *, int, float *);
void   ffft_(float *, int *, int *, int *);
void   fftfit_(float *, float *, float *, int *, float *, float *, float *, 
	       float *, float *, float *, int *);

void   InitPars(struct ASPHdr *);
int    AllocRunMode(struct RunVars *);
int    ReadHdr(struct ASPHdr *, fitsfile *);
int    ReadASPHdr(struct ASPHdr *, fitsfile *);
int    ReadPSRFITSHdr(struct ASPHdr *, fitsfile *);
int    GetTelescope(struct ASPHdr *, struct Telescope *Tel);
void   GetCal(struct ASPHdr *, struct RunVars *, double **);
int    GetCalPhases(double *, int , int *);
double GetCalHeight(double *, int , int *, int *, double *, double *);
int    QuickPlot(double *, int);
int    Median(double *, double *, int, int);
int    ReadData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		fitsfile *, int, int, 
		double **, double **, double **, double **, long **);
int    ReadASPData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		   fitsfile *, int, int, double **, double **, double **, 
		   double **, long **);
int    ReadPSRFITSData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		       fitsfile *, int, int, 
		       double **, double **, double **, double **);
void   MakeStokes(struct ASPHdr *,struct RunVars *, struct StdProfs *,double *,
		  double *, double *, double *, double *);
int    Dedisperse(struct StdProfs *, struct RunVars *, 
		  struct ASPHdr *, struct SubHdr *, 
		  int i_chan);
int    PhaseShift(struct Polyco *, int , struct StdProfs *, struct RunVars *, 
		  struct ASPHdr *, struct SubHdr *, int );
void   BinDown(struct RunVars *, struct StdProfs *, struct StdProfs *);
double DutyLookup(char *);
void   BMask(float *,int *,double *,double *);
void   Baseline(float *,double *, int *, double *, double *);
void   RemoveBase(struct RunVars *, int, struct StdProfs *);
void   MakePol(struct RunVars *, int, struct StdProfs *);
void   ISort(int , int *);
void   FSort(int , float *);
void   DSort(int , double *);
void   FZero(float *, int);
void   DZero(double *, int);
void   IZero(int *, int);
void   LZero(long *, int);
int    ISum(int *, int);
float  FSum(float *, int);
double DSum(double *, int);
int    ArrayZero(float *, int);
double FindPeak(float *, int *, int *);
double Max(double *,int , int *);
double Min(double *,int , int *);
float FMax(float *,int , int *);
float FMin(float *,int , int *);
int    IMax(int *, int, int *);
int    IMin(int *, int, int *);
void   WriteStokes(struct RunVars *, struct StdProfs *, char *, char *);
void   RotateProf(struct RunVars *, struct StdProfs *, double);
void   FitAngle(struct RunVars *, struct ASPHdr *, struct SubHdr *, 
		struct StdProfs *, struct Telescope *);
double GetChi(double, char *, double, double, struct Telescope *);
double ratorad(char *);
double dectorad(char *);
int    WrtASPHdr(struct ASPHdr *, fitsfile *);
int    WrtASPStokes(struct ASPHdr, struct SubHdr, fitsfile *, int, 
		    struct StdProfs *, struct RunVars *);
void   ReadASPStokes(struct ASPHdr *, struct SubHdr *,fitsfile *, long , 
		     struct StdProfs *, int, int);
int    ReadASPAsc(char *, char *, int *, struct StdProfs *, int *);
int    GetThetaBB(struct RunVars *, struct ASPHdr *);
int    FitThetaBB(struct RunVars *, struct ASPHdr *, struct StdProfs *, 
		  int, int);
int    GetMueller(char *, float *, struct ASPHdr *);
int    FitMueller(struct RunVars *, struct ASPHdr *, struct StdProfs *, int);
void   MJDPaste(int, double, char *);
void   TemplateCutoff(struct StdProfs *, int, float);
void   Zero2One(struct StdProfs *, int, char *);
