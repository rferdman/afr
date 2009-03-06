#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


double GetChi(char *Source, double MJD, char *NObs, double RA, double Dec)
{

  int    Verbose=0;
  double Chi;
  double Lat, Long;
  double LST, LSTdeg;
  double HA;

  static int first=1;
  FILE *fptest;


  /* Assign latitudes and longitudes based on telescope -- 
     only at Arecibo and GBT now */

  //  if(NObs==1){
  if(!strcmp(NObs,"1")){
    //    printf("GBT! GBT!\n");fflush(stdout);
    Lat  = 38.0 + 25.0/60.0 + 59.26/3600.0;
    Long = 79.0 + 50.0/60.0 + 23.42/3600.0;
  }
    //  else if(NObs==3){
  else if(!strcmp(NObs,"3")){
    Lat  = 18.0 + 20.0/60.0 + 36.6/3600.0;
    Long = 66.0 + 45.0/60.0 + 11.1/3600.0;  
  }
  else if(!strcmp(NObs,"f")) {  // Nancay...
    Lat = 47.0 + 23.0/60.0;
    Long = -2.0 + 12.0/60.0;   // Nancay is east of GMT
  }
  else{
    printf("The telescope you data comes from does not use ASP or GASP!\n");
    return -999999.;
  }

  /* calculate LST in radians -- based on GST at Dec. 31, 2003 midnight UT */
  LSTdeg = 0.2750375 + 1.0027379093*(MJD-53004.0) - Long/360.0;
  LST = TWOPI*(LSTdeg-floor(LSTdeg));

  if (LST < 0.0)
    LST = TWOPI-LST;

  HA = LST-RA;
  Lat = Lat*TWOPI/360.;

  /*** do we need dec in radians????? ***/

  Chi = 2.0*atan( (sin(HA)*cos(Lat)) / 
		  ( sin(Lat)*cos(Dec) - cos(Lat)*sin(Dec)*cos(HA)) );

  if(Verbose)
    printf("RA = %lf, Dec = %lf, LST = %lf, Chi = %lf\n",RA, Dec, LST, Chi);fflush(stdout);


  /*****  commenting out cretion of pa_test.dat file
    if((fptest = fopen("pa_test.dat","a")) == 0)
      {printf("Can't open PA test file... Exiting...\n");exit(1);}
    fprintf(fptest,"%lf  %lf  %lf  %lf\n",MJD, Chi, LST, HA);fflush(stdout);
    fclose(fptest);
  *****/

  return Chi;



}
