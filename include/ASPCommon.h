/* .h files, define's, and structures used by all programs */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include "ASPHeader.h"
#include "polyco.h"

#define NBINMAX 4096
#define NB      128
#define MAXSPLIT 32
#define NFILEMAX 1000
#define MAXDUMPS 512
#define MAXOMIT  8192
/* #define NCHMAX 32 */

//#define TWOPI 6.2831853071795864
/* MORE DECIMALS...!  */
#define TWOPI 6.2831853071795864769252867665590057683943387987502
#define PI    3.1415926535897932384626433832795028841971693993751

#define DFFAC 2.410e-4


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
  int    FirstChanAdd[NCHMAX];
  int    LastChanAdd[NCHMAX];
  int    MinChans2Add[NCHMAX];
  int    MaxChans2Add[NCHMAX];
  int    NumChans2Add[NCHMAX];
  int    CurZapChan[NCHMAX];
  int    ZapChan[NCHMAX];
  int    DumpOmit[MAXOMIT];
  int    ChanOmit[MAXOMIT];
  int    TotScans[MAXDUMPS*NCHMAX];
  int    TotOmit[MAXDUMPS*NCHMAX];
  int    AllDumpOmit[MAXOMIT];
  int    OmitFlag[MAXDUMPS*NCHMAX];
  float  MM[16*NCHMAX];
  double FSky;
  double ThetaBB[NCHMAX];
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
  float  rstdinv[NBINMAX];
  float  stds[NBINMAX][2];
  float  stdq[NBINMAX][2];
  float  stdu[NBINMAX][2];
  float  stdv[NBINMAX][2];
  float  stdinv[NBINMAX][2];
  float  stdlin[NBINMAX];
  float  stdphi[NBINMAX];
  float  stdphierr[NBINMAX];
  float  stdamp[NBINMAX];
  float  stdpha[NBINMAX];
  float  stdampinv[NBINMAX];
  float  stdphainv[NBINMAX];
  float  stdamps[NBINMAX];
  float  stdampq[NBINMAX];
  float  stdampu[NBINMAX];
  float  stdampv[NBINMAX];
  float  stdphas[NBINMAX];
  float  stdphaq[NBINMAX];
  float  stdphau[NBINMAX];
  float  stdphav[NBINMAX];
  float  profmax[NBINMAX];
  float  profmin[NBINMAX];
  float  stdampx[NBINMAX];
  float  stdampn[NBINMAX];
  float  stdphax[NBINMAX];
  float  stdphan[NBINMAX];
  double SNR;
  double Srms;
  double Qrms;
  double Urms;
  double Vrms;
  double Linrms;
};

/* Functions */
void   cprofc(float *, int, float *, float *);
void   uncprofc(float *, float *, int, float *);
void   ffft_(float *, int *, int *, int *);
void   fftfit_(float *, float *, float *, int *, float *, float *, float *, 
	       float *, float *, float *, int *);

void   InitPars(struct ASPHdr *);
int    ReadASPHdr(struct ASPHdr *, fitsfile *);
void   GetCal(struct ASPHdr *, struct RunVars *, double **);
int    GetCalPhases(double *, int , int *);
double GetCalHeight(double *, int , int *, int *, double *, double *);
int    QuickPlot(double *, int);
int    Median(double *, double *, int, int);
void   ReadASPData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		   fitsfile *, int, long, double **, double **, double **, 
		   double **, int **, char **);
int    ReadStd(struct RunVars *, struct RunVars *, struct StdProfs *, double *);
int    MakeStd(struct ASPHdr *, struct RunVars *, struct StdProfs *, 
	       struct StdProfs *, int);
void   MakeStokes(struct ASPHdr *,struct RunVars *, struct StdProfs *,double *,
		  double *, double *, double *, double *);
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
double FindPeak(float *, int *, int *);
double Max(double *,int , int *);
double Min(double *,int , int *);
int    IMax(int *, int, int *);
int    IMin(int *, int, int *);
void   WriteStokes(struct RunVars *, struct StdProfs *, char *, char *);
void   AddChans(struct RunVars *, struct StdProfs *, struct StdProfs *,
                struct StdProfs *, char **, struct ASPHdr *,
                struct SubHdr *, struct SubHdr *, double *, char **, int);
/*void   AddChans(struct RunVars *, struct StdProfs *, struct StdProfs *, 
		char **, struct ASPHdr *, 
		struct SubHdr *, struct SubHdr *, double *, char **, int);*/
void   RotateProf(struct RunVars *, struct StdProfs *, double);
void   FitAngle(struct RunVars *, struct ASPHdr *, struct SubHdr *, 
		struct StdProfs *);
//double GetChi(char *, double, int, double, double);
double GetChi(char *, double, char *, double, double);
double ratorad(char *);
double dectorad(char *);
void   AddDumps(struct RunVars *, struct ASPHdr *, struct SubHdr *, 
		struct SubHdr *,struct StdProfs *, struct StdProfs *, 
		char **, char **, int, int *);
void   GetAngle(struct RunVars *, struct StdProfs *, struct StdProfs *);
void   FinalDump(struct RunVars *, struct ASPHdr *, struct StdProfs *, 
		 char **);
int    WrtASPHdr(struct ASPHdr *, fitsfile *);
/* int    WrtASPStokes(struct ASPHdr, struct SubHdr, fitsfile *, int, 
   struct StdProfs *, Cmdline *, struct RunVars *); */
int    WrtASPStokes(struct ASPHdr, struct SubHdr, fitsfile *, int, 
		    struct StdProfs *,  int, struct RunVars *);
void   ReadASPStokes(struct ASPHdr *, struct SubHdr *,fitsfile *, long , 
		     struct StdProfs *, int, int);
int    ReadASPAsc(char *, char *, int *, struct StdProfs *, int *);
int    GetThetaBB(struct RunVars *, struct ASPHdr *);
int    FitThetaBB(struct RunVars *, struct ASPHdr *, struct StdProfs *, 
		  int, int);
int    GetMueller(char *, float *, struct ASPHdr *);
int    FitMueller(struct RunVars *, struct ASPHdr *, struct StdProfs *, int);
void    MJDPaste(int, double, char *);
