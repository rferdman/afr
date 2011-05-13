/* Calculate difference between current pulse phase and expected phase, if
   using a new and improved polyco.dat file to do so.  Useful for re-aligning
   profiles that have been drifting over time due to older polyco files */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"
#include "phase_calc.h"

int PhaseShift(struct Polyco *Polycos, int n_poly, struct StdProfs *Profile, 
	       struct RunVars *RunMode, struct ASPHdr *Hdr, 
	       struct SubHdr *SubHdr, int chan){

  
  int    poly_set_used;
  double DumpMiddleDays, RefPhase, RefFreq;
  double PhaseDiff, ByAngle;

  DumpMiddleDays = SubHdr->DumpMiddleSecs /86400.;

  if((poly_set_used = PhaseCalc(Polycos, n_poly, Hdr->obs.IMJDStart, 
				DumpMiddleDays, &RefPhase, &RefFreq)) < 0){
    printf("Could not find any polyco sets to match MJD. Exiting...\n");
    return -1;
  }
  
  
  PhaseDiff = RefPhase - SubHdr->DumpRefPhase[chan];
  ByAngle = PhaseDiff*TWOPI;
  //  ByAngle = PhaseDiff;
    
  RotateProf(RunMode, Profile, ByAngle);

  printf("%.5lf : oldphase=%.5lf  newphase=%.5lf\n", Hdr->obs.ChanFreq[chan], SubHdr->DumpRefPhase[chan], RefPhase);
  
  SubHdr->DumpRefPeriod[chan] = 1./RefFreq;
  SubHdr->DumpRefPhase[chan]  = RefPhase;

  return poly_set_used;


}
