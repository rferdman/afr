
#ifndef _POLYCO_H
#define _POLYCO_H

#define MAX_PC_COEFF 15
#define MAX_PC_SETS  2048

//#include "asp_params.h"
// #include "misc.h"

struct Polyco {
  int NSets;                   /* Number of sets for this pulsar */
  int NMinutes;                /* Number of minutes' validity */
  int NCoeff;                  /* Number of coefficients (max 15) */
  double DM;                   /* DM of this pulsar */
  double FSkyRef;              /* Reference rotational frequency */
  double FRotRef;              /* Reference rotational frequency */
  double PhRotRef;             /* Reference rotational phase */
  double EarthZ4;              /* Earth motion correction to DM */
  double MjdMidInt;            /* Integer part of MJD for this set */
  double MjdMidFrac;           /* Fractional part of MJD for this set */
  double Coeff[MAX_PC_COEFF]; /* Actual coefficients */
};

int GetPoly(char *, char *, struct Polyco *, double, double);
//int GetPoly(struct asp_params *obs_params, struct Polyco *pc, double ChanFreq);
#endif

