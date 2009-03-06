#include <stdio.h>
#include <math.h>
#include "cpgplot.h"
#include "fitsio.h"
#include "ASPCommon.h"

/* Various tools for retrieving and using calibration data */


/* Find phases at which cal pulse is on/off */
int QuickPlot(double *Profile, int NPtsProf)
{
  int   i, retval=0;
  int   iymin, iymax;
  float xmin,  xmax,  ymin,  ymax;
  float Prof[NBINMAX], Bins[NBINMAX];

  cpgbeg(0, "/xs", 1, 1);


  /* Get min and max values of x and y axes */
  xmin = 0.;
  xmax = (float)NPtsProf;
  ymin = (float)(Min(Profile, NPtsProf, &iymin));
  ymax = (float)(Max(Profile, NPtsProf, &iymax));

  cpgenv(xmin, xmax, ymin, ymax, 0, 0);
  cpglab("Bin", "Counts", "Profile Plot");

  for(i=0;i<NPtsProf;i++){
    Prof[i] = (float)Profile[i];
    Bins[i] = (float)i;
  }

  cpgline(NPtsProf, Bins, Prof);
  cpgend();

  return retval;
}
