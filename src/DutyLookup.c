/* Looks up duty cycle for different pulsars */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ASPCommon.h"

double DutyLookup(char *Source)
{

  static int first=1;
  int i_pt = 0;
  double Duty;

  Duty = 0.95;     /*  default  */


  if(strncmp(&Source[i_pt],"c",1) == 0) Duty = 0.45;  /* Cal scan case */

  /* If there is a "B" or "J" at the front of the name, move name pointer up by one */
  if(strncmp(&Source[0],"J",1)==0 || strncmp(&Source[0],"B",1)==0) i_pt++;

  if(strncmp(&Source[i_pt],"0034-05",7) == 0) Duty = 0.50;
  if(strncmp(&Source[i_pt],"0218+42",7) == 0) Duty = 0.30;
  if(strncmp(&Source[i_pt],"0613-",5) == 0) Duty = 0.63;
  if(strncmp(&Source[i_pt],"0621+10",7) == 0) Duty = 0.60;
  if(strncmp(&Source[i_pt],"0737-3039A",10) == 0) Duty = 0.53; //estimation
  if(strncmp(&Source[i_pt],"0751+18",7) == 0) Duty = 0.70;
  if(strncmp(&Source[i_pt],"1012+53",7) == 0) Duty = 0.40;
  if(strncmp(&Source[i_pt],"1022+10",7) == 0) Duty = 0.83;
  if(strncmp(&Source[i_pt],"1023+0038",9) == 0) Duty = 0.50; // estimation IHS
  if(strncmp(&Source[i_pt],"1257+12",7) == 0) Duty = 0.65;
  if(strncmp(&Source[i_pt],"1518+49",7) == 0) Duty = 0.85;
  if(strncmp(&Source[i_pt],"1534+12",7) == 0) Duty = 0.40;
  if(strncmp(&Source[i_pt],"1620-26",7) == 0) Duty = 0.80;
  if(strncmp(&Source[i_pt],"1643-12",7) == 0) Duty = 0.75;
  if(strncmp(&Source[i_pt],"1713+07",7) == 0) Duty = 0.70;
  if(strncmp(&Source[i_pt],"1730-23",7) == 0) Duty = 0.78;
  if(strncmp(&Source[i_pt],"1744-11",7) == 0) Duty = 0.93;
  if(strncmp(&Source[i_pt],"1744-24",7) == 0) Duty = 0.91;
  if(strncmp(&Source[i_pt],"1756-2251",9) == 0) Duty = 0.80; // estimation RDF
  if(strncmp(&Source[i_pt],"1820-30",7) == 0) Duty = 0.83;
  if(strncmp(&Source[i_pt],"1802-2124",9) == 0) Duty = 0.90; // estimation RDF
  if(strncmp(&Source[i_pt],"1821-24",7) == 0) Duty = 0.50;
  if(strncmp(&Source[i_pt],"1855+09",7) == 0) Duty = 0.50;
  if(strncmp(&Source[i_pt],"1911-11",7) == 0) Duty = 0.91;
  if(strncmp(&Source[i_pt],"1913+16",7) == 0) Duty = 0.70;
  if(strncmp(&Source[i_pt],"1929+10",7) == 0) Duty = 0.90;
  if(strncmp(&Source[i_pt],"1937+21",7) == 0) Duty = 0.91;
  if(strncmp(&Source[i_pt],"1957+20",7) == 0) Duty = 0.65;
  if(strncmp(&Source[i_pt],"2019+24",7) == 0) Duty = 0.70;
  if(strncmp(&Source[i_pt],"2033+17",7) == 0) Duty = 0.60;
  if(strncmp(&Source[i_pt],"2051-08",7) == 0) Duty = 0.78;
  if(strncmp(&Source[i_pt],"2145-07",7) == 0) Duty = 0.70;
  if(strncmp(&Source[i_pt],"2229+26",7) == 0) Duty = 0.55;
  if(strncmp(&Source[i_pt],"2317+14",7) == 0) Duty = 0.55;

  
  if(first){
    printf("Duty cycle for %s = %lf\n\n",Source,Duty);fflush(stdout);
    first=0;
  }

  return(Duty);

}
