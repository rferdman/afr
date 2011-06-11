#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

double GetLST(double MJD, double Longitude)
{
  int    Verbose=0;
  double longi,lstrad;
  double a,b,c,d;
  double seconds_per_jc, s2r;
  double bprime,cprime,dprime;
  double nmjdut,fmjdut,mjdut;
  double tu0,dtu,tu,gmst0,sdd,gst,xlst; //*test1,*test2,*test3,*test4;


  /* If I don't allocate memory for all these variables other parts of the code seem to not like it very much */
/*   lstrad = (double *)malloc(sizeof(double));
  longi = (double *)malloc(sizeof(double));
  nmjdut = (double *)malloc(sizeof(double));
  fmjdut = (double *)malloc(sizeof(double));
  mjdut = (double *)malloc(2*sizeof(double));
  tu0 = (double *)malloc(sizeof(double));
  dtu = (double *)malloc(sizeof(double));
  tu = (double *)malloc(sizeof(double));
  gmst0 = (double *)malloc(sizeof(double));
  sdd = (double *)malloc(sizeof(double));
  gst = (double *)malloc(sizeof(double));
  xlst = (double *)malloc(sizeof(double)); */
  /** test1 = (double *)malloc(sizeof(double));
  test2 = (double *)malloc(sizeof(double)); 
  test3 = (double *)malloc(sizeof(double)); 
  test4 = (double *)malloc(sizeof(double)); **/

  a = 24110.54841;
  b = 8640184.812866;
  c = 0.093104;
  d = -6.2e-6;
  s2r = 7.272205216643039903848711535369e-5;

  seconds_per_jc = 86400.*36525;

  bprime = 1. + b/seconds_per_jc;
  cprime = 2. * c/seconds_per_jc;
  dprime = 3. * d/seconds_per_jc;  

  nmjdut=floor(MJD); //The integer part of the MJD
  fmjdut=MJD-nmjdut; //The fractional part of the MJD

  /*if (Verbose)
    printf("The old MJD is %lf %lf \n",*nmjdut,*fmjdut);fflush(stdout); 

  Adjust MDJ for leap seconds
  GetLeaps(*nmjdut,*fmjdut,mjdut);

  *nmjdut=mjdut[0]; 
  *fmjdut=mjdut[1];*/
  

  if (Verbose)
    printf("The MJD is %lf %lf \n",nmjdut,fmjdut);fflush(stdout);

  
  /**tu0 = ((*nmjdut-51545.)+0.5)/3.6525e4;*/
  tu0 = (nmjdut-51544.5)/3.6525e4;
  dtu = fmjdut/3.6525e4;
  tu = tu0 + dtu;


  /* This is from TEMPO */ 
  gmst0 = (a + tu0*(b + tu0*(c + tu0*d)))/86400.;

  sdd = bprime + tu*(cprime + tu*dprime);
  
/* Next line was incorrect -- June 8, 2011 */
/*   gst=gmst0 + (seconds_per_jc/86400. + b/86400. + d*(tu*tu + tu*tu0 + tu0*tu0)/86400.)*dtu; */
  gst=gmst0+dtu*(seconds_per_jc + b + c*(tu+tu0) + d*(tu*tu + tu*tu0 + tu0*tu0))/86400.;
 

  /* This is from SLALIB */
  /* *gst = ((a + *tu*(b + *tu*(c + *tu*d)))/86400) + *dtu; */

  /*** This gives LST in number of turns ***/
  longi=Longitude/360.;
  xlst=(gst-longi)-floor(gst-longi);

  /**xlst=xlst-floor(xlst); **/

  if(Verbose)
    printf("tu0 = %e, dtu = %e, tu = %e, gmst0 = %e, gst = %e, longi = %lf, xlst = %lf Longitude = %f\n",tu0,dtu,tu,gmst0,gst,longi,xlst, Longitude);fflush(stdout);

  /*** Change to Radians ***/
  lstrad=TWOPI*xlst;

  if(lstrad < 0.){
    lstrad = TWOPI-lstrad;
  }

  /*lstrad=*lst;
  LST = lstrad;*/

  if(Verbose)
    printf("LST in radians is %lf \n",lstrad);fflush(stdout);

/*  free(lstrad); 
  free(longi); 
  free(nmjdut); 
  free(fmjdut); 
  free(tu0); 
  free(dtu); 
  free(tu);
  free(gmst0);
  free(sdd);
  free(gst);
  free(xlst); */
  /** free(test1);
  free(test2);
  free(test3);
  free(test4); **/
  
  return lstrad;
  
}
