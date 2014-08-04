
/* ========================================================================= */
/* Header parameter file for ASP machine                                     */
/*            R. Ramachandran, 22/Jan/03, Berkeley, CA.                      */
/* ========================================================================= */
#include "ASPDefs.h"
#include "polyco.h"


struct General
{ char     ScanName[STRINGLEN];   /* Scan Name */
  char     SoftVer[STRINGLEN];    /* Software version number */
  char     FitsType[STRINGLEN];   /* FITS Descrition type */
  char     Platform[STRINGLEN];   /* Platform of the software */
  char     HdrVer[STRINGLEN];     /* Header version */
  char     Observer[STRINGLEN];   /* Observer's name 15 characters! */
  char     ProjID[STRINGLEN];     /* Project ID */
  char     CommentOper[120];      /* Comments by Operator during Observation */
  char     FEName[STRINGLEN];     /* Front end name */
  char     FEPol[STRINGLEN];      /* Front end polarisation nature LIN/CIRC */
  char     IBeam[STRINGLEN];      /* Beam ID for multiple-beam systems */
  char     PntID[STRINGLEN];      /* Name or ID for pointing centre (multibeam feeds) */
  char     BEName[STRINGLEN];     /* Backend name */
  char     BEConfFile[STRINGLEN]; /* Backend config file */
  char     ObsMode[STRINGLEN];    /* Observation Mode PSR/CAL/STRCTRAL */
  char     FeedMode[STRINGLEN];   /* Feed track mode -- FA/CPA/SPA/TPA */
  char     CalMode[STRINGLEN];    /* Cal mode (OFF, SYNC, EXT1, EXT2) */
  int      NRcvr;                 /* Number of receiver polarization channels */
  int      FDHand;                /* Hand of polarization: +/- 1.  +1 is LIN: A=X, B=Y, CIRC:A=L, B=R */
  float    FDSang;                /* FA of E vector for equal sig in A&B */
  float    FDXYph;                /* Phase of A*B for injected cal */ 
  double   FAReq;                 /* Feed/position angle requested */
  int      BEPhase;               /* BE cross-phase: 0 unknown, +/-1 std/rev */
  int      BEDownConv;            /* BE downconversion conjugation */     
  float    BEDelay;               /* Backend propn delay from digitiser in seconds */
  double   TCycle;                /* On-line cycle time */
};

struct Target
{
  char     PSRName[STRINGLEN];    /* Source Name */
  double   RA;                    /* Right Assention hhmmss.ssss */
  double   Dec;                   /* Declination ddmmss.ssss */

  double   BeamMajor;             /* Beam major axis length */
  double   BeamMinor;             /* Beam minor axis length */
  double   BeamPA;                /* Beam position angle */

  double   Ant_X;                 /* Antenna ITRF X-coordinate in m */
  double   Ant_Y;                 /* Antenna ITRF Y-coordinate in m */
  double   Ant_Z;                 /* Antenna ITRF Z-coordinate in m */

  float    Epoch;                 /* Coordinate Epoch 1950.0/2000.0 */
  char     CoordMode[STRINGLEN];  /* Coordinate mode for the following 
				     4 entries */
  char     StartCrd1[STRINGLEN];  /* Coordinate 1 at start time */
  char     StartCrd2[STRINGLEN];  /* Coordinate 1 at start time */
  char     StopCrd1[STRINGLEN];   /* Coordinate 1 at stop time */
  char     StopCrd2[STRINGLEN];   /* Coordinate 1 at stop time */
  double   StartLST;              /* LST at the start of observation */
  char     TrackMode[STRINGLEN];  /* Track mode TRACK/SLEW/DRIFT etc. */
};

struct ObservationMode
{
  int      IMJDStart;             /* Start Day number (MJD) */
  int      StartTime;             /* Start time (since 0.0 hr, in sec) */
  double   StartFSec;             /* Fractional sec after StartTime */
  char     StartDate[STRINGLEN];  /* Start Date */
  char     StartUT[STRINGLEN];    /* Start Time UT */
  int      NSubOffs;              /* Number of subintegrations current data file is offset from beginning of observation */
  double   ObsLength;             /* Length of Observation (sec) */
  float    IonRM;                 /* Ionospheric Rotation measure correction */
  double   ClockOffset;           /* [Obs.Std - GPS] (microsec) */
  char     ObsvtyCode[16];         /* Observatory code (TEMPO-compatible) */
  int      NChanOrig;             /* Original no. of frequency channels */
  int      NChan;                 /* No. of frequency channels */
  int      NPoln;                 /* No. of polarizations */
  char     PolnType[STRINGLEN];   /* Polarization identifier (eg. AABBCRCI, IQUV, etc.) */
  double   FSkyCent;              /* Overall sky centre frequency */
  double   ChanFreqMean;          /* Calculated mean of min and max channel frequencies */
  double   BW;                    /* Observing bandwidth */
  int      ChanIndx[NCHMAX];      /* Orig channel index in xASP */
  double   ChanFreq[NCHMAX];      /* Sky mid-frequency (MHz) */
  double   ChanWidth;             /* Frequency Channel width (MHz) */
  float    ChanTSys[NCHMAX];      /* System Temperature (K) of channels */
  float    ChanTCal[NCHMAX];      /* Calibration temperature (K) of channels */
  double   CalFreq;               /* Cal modulation frequency in Hz */
  double   CalDutyCycle;          /* Cal duty cycle */
  double   CalPhase;              /* Cal phase (wrt start time) */
  double   CalNStates;            /* Number of states in cal pulse */
  int      Sideband;              /* Forward/Reverse band? 
				     [1 -> forward / -1 -> reverse] */
  double   DM;                    /* Dispersion measure (pc/cc) */
  double   ChanDM;                /* Dispersion measure used for online dedispersion (pc/cc) */
  int      DMMethod;              /* [0,no DDP], [1,CohDD], [2,IncohDD] */
  int      ChanFFTLen[NCHMAX];    /* log_2[NFFT] for each channel for Chh.DD */
  int      ChanChirpLen[NCHMAX];  /* Chirp length in each channel */
  int      ChanOverlap[NCHMAX];   /* Amount of convolution overlap 
				     in each channel */
  double   RM;                    /* Rotation measure */
  int      RMMethod;              /* [0,not applied], [1,coherent], 
				     [2,incoherent] */
  double   SampInterval;          /* sampling interval (seconds) */
  char     OPString[10];          /* Output string, I.Q.U.V.L.T, etc. */
  int      OPScale;               /* Output scale: [0,no cal], [1,K], [2,Jy] */
  Boolean  DynDC;                 /* Mean/RMS Dynamically corrected? 
				     [TRUE/FALSE] */
  float    DynTime;               /* Time scale assumed for Dyn. 
				     mean/rms correction (sec) */
  int      NBitPerSamp;           /* No. of bits per sample */
  Boolean  InterfString;          /* Interference corrected?  [TRUE/FALSE] */
  char     ZapFileName[40];       /* Interference Zap filename */
};

struct ASPSetup
{
  int      NDS;                   /* No. of data servers */
  char     ConfFile[STRINGLEN];   /* ABPP Config file [like mf_conf.dat] */
  float    IFLO[NDSMAX];          /* Intermediate frequency LO values */
  float    RFLO[NDSMAX];          /* Radio Frequency LO values */
  float    MFAttn[2*NDSMAX];      /* MF Attenuation factors [X,Y] */
  float    MFFilter;              /* Low pass filter value (Hz) */
  float    MFPower[2*NDSMAX];     /* Power out of MF [volts] */
  double   DFBRate;               /* A/D sampling rate in DFB */
  float    DFBLO[NCHMAX];         /* Digital mixing by LOs in DFB */
  float    DFBGain[NCHMAX];       /* Gain factor, Harris' FIR coeffecients */
  int      DFBDecimate;           /* DFB Decimation factor (FIR factor) */
  int      DFBMode;               /* EPROM mode */
  int      FIRFactor;             /* FIR-factor, if any */
  double   SampleRate;            /* [DFBRate/DFBDecimate] (Hz) */
  int      SplitFac;              /* ABPP channel split factor */
};

struct FileContent
{
  int      SizeHdr;                   /* Size of header block (bytes) */
  int      SizeDataBlk;               /* Size of data block (bytes) */
  int      IMJDFirstSamp;             /* Integer day number of first sample (MJD) */
  double   fMJDFirstSamp;             /* fractional day (since 0.0 hrs) 
					 for first sample */
  char     FileDate[STRINGLEN];  /* File creation date (YYYY-MM-DD) */
  char     FileUT[STRINGLEN];    /* File creation time (UT) */ 
};

struct Reduction
{
  int      RIMJDFirstSamp;        /* Day number of First sample (MJD) */
  double   RfMJDFirstSamp;        /* fractional day since 0.0 hrs for 
				     first sample */
  int      RNTimeDumps;           /* No. of time dumps 
				     [timedump = folded prof; t-series etc.] */
  double   TDump;                 /* Length of each time dump */
  int      RNBinTimeDump;         /* No. of bins in each time dump */
  int      RChanNum;              /* No. of frequency channels */
  double   RChanWidth;            /* Width of each frequency channel (MHz) */
  char     ROPString[STRINGLEN];  /* Output string I.Q.U.V.L.T, etc. */
  double   RDM;                   /* Assumed DM (pc/cc) */
  float    RRM;                   /* Assumed RM */
  double   RPeriod;               /* Rotation period at first sample (sec) */
  char     RFoldType;             /* Type of folding 
				     [1,polyco], [2,ConstPeriod] */
  int      RNPolyCoeff;           /* No. of Polyco coeffecients */
  char     RDateString[40];       /* Date and time of analysis */
  char     RUserName[16];         /* Analyser's name */
  int      ROPModeCode;           /* Output mode code (binary-coded) */
  char     RCommString[120];      /* Command string given for analysis */
/* Data has been dedispersed (HISTORY table in PSRFITS files) */
  int      Dedisp;                 
  /* Include polyco here for PSRFITS data, from which we read (or create) polycos
     corresponding to original observation, of use when realigning profiles.  Will 
     allocate later, only once we know that we need them */
  int      NPoly;   /* Number of polyco sets found in header or from polycos generated from parfile */
  struct Polyco *Polycos; /* Be sure to allocate the memory for this... */
};

struct ASPHdr
{
  struct General          gen;
  struct Target           target;
  struct ObservationMode  obs;
  struct ASPSetup         asp;
  struct FileContent      filecont;
  struct Reduction        redn;
};


struct SubHdr
{
  double DumpMiddleSecs; 
  double DumpRefPeriod[NCHMAX];
  double DumpRefPhase[NCHMAX];
};
