/* Fills in info about telescope corresponding to current data file

   -- R. Ferdman, 7 May 2010 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"


int GetTelescope(struct ASPHdr *hdr, struct Telescope *Tel)
{

  /* This seems to get rid of leading or trailing white spaces... */
  sscanf(hdr->obs.ObsvtyCode, "%s", hdr->obs.ObsvtyCode);

  /* GBT */
  if(!strcmp(hdr->obs.ObsvtyCode, "1") 
     || !strncasecmp(hdr->obs.ObsvtyCode, "GBT", 3) ){
    sprintf(hdr->obs.ObsvtyCode, "1");
    Tel->Lat  = 38.0 + 25.0/60.0 + 59.26/3600.0;
    Tel->Long = 79.0 + 50.0/60.0 + 23.42/3600.0;
  }
  /* Arecibo */
  else if(!strcmp(hdr->obs.ObsvtyCode, "3") 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "ARECIBO", 7)
 	  || !strncasecmp(hdr->obs.ObsvtyCode, "AO", 2)){
   sprintf(hdr->obs.ObsvtyCode, "3");
    Tel->Lat  = 18.0 + 20.0/60.0 + 36.6/3600.0;
    Tel->Long = 66.0 + 45.0/60.0 + 11.1/3600.0;  
  }
  /* Nancay */
  else if(!strcmp(hdr->obs.ObsvtyCode, "f") 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "NANCAY", 6) 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "NCY", 3)) {  
    sprintf(hdr->obs.ObsvtyCode, "f");
    Tel->Lat  = 47.0 + 23.0/60.0;
    Tel->Long = -(2.0 + 12.0/60.0);   // Nancay is east of GMT
  }
  /* Parkes */
  else if(!strcmp(hdr->obs.ObsvtyCode, "7") 
	  || !strcasecmp(hdr->obs.ObsvtyCode, "PARKES") 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "PKS", 3)) {  
    sprintf(hdr->obs.ObsvtyCode, "7");
    Tel->Lat  = -32.0 + 59.0/60.0  + 54.263/3600.0;
    Tel->Long = 148.0 + 15.0/60.0  + 48.636/3600.0;  
  }
  /* Jodrell Lovell */
  else if(!strcmp(hdr->obs.ObsvtyCode, "8") 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "JODRELL", 7) 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "JB", 3)  
	  || !strncasecmp(hdr->obs.ObsvtyCode, "LOV", 3)) {  
    sprintf(hdr->obs.ObsvtyCode, "8");
    Tel->Lat  = 53.0 + 14.0/60.0 + 10.50/3600.0 ;
    Tel->Long =  2.0 + 18.0/60.0 + 25.74/3600.0 ;  
  }
  /* Effelsberg */
  else if(!strcmp(hdr->obs.ObsvtyCode, "g") 
	  || !strncasecmp(hdr->obs.ObsvtyCode, "EFF", 3)) {  
    sprintf(hdr->obs.ObsvtyCode, "8");
    Tel->Lat  = 50.0 + 31.0/60.0 + 30.0/3600.0 ;
    Tel->Long = -(6.0 + 53.0/60.0 + 0.3/3600.0);  // Effelsberg is east of GMT
  }
  else{
    printf("The telescope your data comes from (%s) is not compatible ",
	   hdr->obs.ObsvtyCode);
    printf("with AFR...\n");
    return -999999.;
  }





  return 1;
}
