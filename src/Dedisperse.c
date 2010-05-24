/* Dedisperse input profile to centre frequency.

   -- R. Ferdman, 13 May 2010  */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

int Dedisperse (struct StdProfs *Profile, struct RunVars *RunMode, 
		struct ASPHdr *Hdr, struct SubHdr *SubHdr, 
		int i_chan)
{


  
  double  dm_delay, phase_delay, ByAngle;


  /* Calculate dispersion delay realtive to centre frequency */
  dm_delay = (Hdr->obs.DM / DFFAC) * 
    ( (1.0/(Hdr->obs.ChanFreq[i_chan]*Hdr->obs.ChanFreq[i_chan])) - 
      (1.0/(Hdr->obs.FSkyCent*Hdr->obs.FSkyCent)) );
  phase_delay = (dm_delay/SubHdr->DumpRefPeriod[i_chan]);
  SubHdr->DumpRefPhase[i_chan] -= phase_delay;
  ///  RotateProf
  //      printf("Chan%d:\nPeriod =  %.4lf\nPhase = %.4lf\ndm_delay = %.2lf\nphase_delay = %.4lf\n",i_chan,SubInHdr[i_dump].DumpRefPeriod[i_chan], SubInHdr[i_dump].DumpRefPhase[i_chan], dm_delay, phase_delay);

  /* Now rotate prof according to dispersion, with reference to 
     centre frequency */
  RunMode->NBins = Hdr->redn.RNBinTimeDump;
  ByAngle = -phase_delay*TWOPI;
  RotateProf(RunMode, Profile, ByAngle);










      return 0;


}