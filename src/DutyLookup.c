/* Looks up duty cycle for different pulsars */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ASPCommon.h"

double DutyLookup(char *Source)
{

  static int first=1;
  double Duty;

  Duty = 0.95;     /*  default  */


  if(strncmp(Source,"c",1) == 0) Duty = 0.45;  /* Cal scan case */

  if(strncmp(Source,"0034-05",7) == 0) Duty = 0.50;
  if(strncmp(Source,"0218+42",7) == 0) Duty = 0.30;
  if(strncmp(Source,"0613-",5) == 0) Duty = 0.63;
  if(strncmp(Source,"0621+10",7) == 0) Duty = 0.60;
  if(strncmp(Source,"0737-3039A",10) == 0) Duty = 0.53; //estimation
  if(strncmp(Source,"0751+18",7) == 0) Duty = 0.70;
  if(strncmp(Source,"1012+53",7) == 0) Duty = 0.40;
  if(strncmp(Source,"1022+10",7) == 0) Duty = 0.83;
  if(strncmp(Source,"1023+0038",9) == 0) Duty = 0.50; // estimation IHS
  if(strncmp(Source,"1257+12",7) == 0) Duty = 0.65;
  if(strncmp(Source,"1518+49",7) == 0) Duty = 0.85;
  if(strncmp(Source,"1534+12",7) == 0) Duty = 0.40;
  if(strncmp(Source,"1620-26",7) == 0) Duty = 0.80;
  if(strncmp(Source,"1643-12",7) == 0) Duty = 0.75;
  if(strncmp(Source,"1713+07",7) == 0) Duty = 0.70;
  if(strncmp(Source,"1730-23",7) == 0) Duty = 0.78;
  if(strncmp(Source,"1744-11",7) == 0) Duty = 0.93;
  if(strncmp(Source,"1744-24",7) == 0) Duty = 0.91;
  if(strncmp(Source,"1756-2251",9) == 0) Duty = 0.80; // estimation RDF
  if(strncmp(Source,"1820-30",7) == 0) Duty = 0.83;
  if(strncmp(Source,"1802-2124",9) == 0) Duty = 0.90; // estimation RDF
  if(strncmp(Source,"1821-24",7) == 0) Duty = 0.50;
  if(strncmp(Source,"1855+09",7) == 0) Duty = 0.50;
  if(strncmp(Source,"1911-11",7) == 0) Duty = 0.91;
  if(strncmp(Source,"1913+16",7) == 0) Duty = 0.70;
  if(strncmp(Source,"1929+10",7) == 0) Duty = 0.90;
  if(strncmp(Source,"1937+21",7) == 0) Duty = 0.91;
  if(strncmp(Source,"1957+20",7) == 0) Duty = 0.65;
  if(strncmp(Source,"2019+24",7) == 0) Duty = 0.70;
  if(strncmp(Source,"2033+17",7) == 0) Duty = 0.60;
  if(strncmp(Source,"2051-08",7) == 0) Duty = 0.78;
  if(strncmp(Source,"2145-07",7) == 0) Duty = 0.70;
  if(strncmp(Source,"2229+26",7) == 0) Duty = 0.55;
  if(strncmp(Source,"2317+14",7) == 0) Duty = 0.55;

  
  if(first){
    printf("Duty cycle for %s = %lf\n\n",Source,Duty);fflush(stdout);
    first=0;
  }

  return(Duty);

}
