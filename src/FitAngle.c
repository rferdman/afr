#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ASPCommon.h"

double ratorad(char *ra)
{
 
  int    hh,mm;
  double ss,rads;
 
  sscanf(ra,"%2d:%2d:%lf",&hh,&mm,&ss);
  rads = (3600.*hh + 60.*mm + ss)*TWOPI/86400.;
  return(rads);
}

double dectorad(char *dec)
{
 
  int dd,mm;
  double ss,rads;
 
  sscanf(dec,"%2d:%2d:%lf",&dd,&mm,&ss);
  rads = (3600.*dd + 60.*mm + ss)*TWOPI/1296000.;
  return(rads);
}

/* Rotate profile using parallactic angle corrections */

void FitAngle(struct RunVars *RunMode, 
	      struct ASPHdr *hdr, struct SubHdr *subhdr,
	      struct StdProfs *Profile,
	      struct Telescope *Tel)
{

  int i;//, NObs;
  double MJD, RA, Dec, AngRot;
  double Wgt;
  struct StdProfs TempProf;
  
  char NObs[3];
  
  if(hdr->target.RA > 0. && fabs(hdr->target.Dec) > 0. ){
    RA = (hdr->target.RA/24.)*TWOPI;
    Dec = hdr->target.Dec*TWOPI/360.;
  }
  else{
    /* Quick hardcoded lookup table for RA */

    if(!strncmp(RunMode->Source,"0531",4)){
      RA  = ratorad("05:34:31.973");
      Dec = dectorad("22:00:52.06");
    }  else if(!strncmp(RunMode->Source,"0751",4)){
      RA  = ratorad("07:51:09");
      Dec = dectorad("18:07:39");
    }  else if(!strncmp(RunMode->Source,"0950",4)){
      RA  = ratorad("09:53:09.309");
      Dec = dectorad("07:55:35.75");
    }  else if(!strncmp(RunMode->Source,"1012+53",7)){
      RA  = ratorad("10:12:33.43370");
      Dec = dectorad("53:07:02.5884");
    }  else if(!strncmp(RunMode->Source,"1133",4)){
      RA  = ratorad("11:36:03.22");
      Dec = dectorad("15:51:05.72");
    }  else if(!strncmp(RunMode->Source,"1534",4)){
      RA  = ratorad("15:37:10");
      Dec = dectorad("11:55:56");
    }  else if(!strncmp(RunMode->Source,"1853+1308",9)){
      RA  = ratorad("18:53:56");
      Dec = dectorad("13:04:53");
    }  else if(!strncmp(RunMode->Source,"1730-2304",9)){
      RA  = ratorad("17:30:21.6483");
      Dec = dectorad("-23:04:31.4");
    } 
    else if(!strncmp(RunMode->Source,"1855",4)){
      RA  = ratorad("18:57:36");
      Dec = dectorad("09:43:17");
    }
    else if(!strncmp(RunMode->Source,"1904+0412",9)){
      RA  = ratorad("19:04:31.39");
      Dec = dectorad("04:12:05.77");
    }
    else if(!strncmp(RunMode->Source,"1905+0400",9)){
      RA  = ratorad("19:05:28.27");
      Dec = dectorad("04:00:10.93");
    }
    else if(!strncmp(RunMode->Source,"1910+1257",9)){
      RA  = ratorad("19:10:09.70");
      Dec = dectorad("12:56:25.53");
    }
    else if(!strncmp(RunMode->Source,"1937",4)){
      RA  = ratorad("19:39:38.56");
      Dec = dectorad("21:34:59.14");
    }
    else{
      printf("Could not find RA for PSR %s! Exiting...\n",RunMode->Source);
      exit(1);
    }

  }
  MJD = (double)hdr->obs.IMJDStart + (subhdr->DumpMiddleSecs)/86400. ;
/*   MJD = (double)hdr->obs.IMJDStart + ((double)hdr->obs.StartTime)/86400.; */

  if(RunMode->Verbose)
    printf("MJDStart = %d, DumpMiddleSecs = %lf\n",hdr->obs.IMJDStart, 
	   subhdr->DumpMiddleSecs); 


  /*  printf("RAIn = %lf, RA = %lf, DECIn = %lf, DEC = %lf\n",hdr->target.RA,RA,
      hdr->target.Dec,Dec);fflush(stdout); */

  //  sscanf(hdr->obs.ObsvtyCode,"%d",&NObs);
  // sscanf(hdr->obs.ObsvtyCode,"%s",NObs);

  // printf("NObs = %d\n",NObs);fflush(stdout);


  AngRot = -1.0*GetChi(MJD, hdr->obs.ObsvtyCode, RA, Dec, Tel);
  //  getchi_(RunMode->Source,&MJD,&ChiInt,&TSky,&NObs,&RA,&Dec,8L);

  /* Assume Circular feeds for now */

/*   AngRot = -ChiInt; */
//   AngRot = ChiInt; 

  if(RunMode->Verbose)
    printf("AngRot = %lf radians = %lf degrees\n",AngRot,AngRot*360./TWOPI);

  /* Assume equal weight for now */

  Wgt = 1.0;

  for(i=0;i<RunMode->NBins;i++){
    TempProf.rstds[i] = Wgt*Profile->rstds[i];
    TempProf.rstdq[i] = Wgt*(cos(AngRot)*Profile->rstdq[i] +
			     sin(AngRot)*Profile->rstdu[i]);
    TempProf.rstdu[i] = Wgt*(cos(AngRot)*Profile->rstdu[i] -
			     sin(AngRot)*Profile->rstdq[i]);
    TempProf.rstdv[i] = Profile->rstdv[i];   
  }

#if 0
  for(i=0;i<RunMode->NBins;i++){
    TempProf.rstds[i] = Wgt*Profile->rstds[i];
    TempProf.rstdu[i] = Wgt*(cos(AngRot)*Profile->rstdu[i] +
			     sin(AngRot)*Profile->rstdq[i]);
    TempProf.rstdq[i] = Wgt*(cos(AngRot)*Profile->rstdq[i] -
			     sin(AngRot)*Profile->rstdu[i]);
    
    TempProf.rstdv[i] = Profile->rstdv[i];   
  }
#endif  

#if 0
  for(i=0;i<RunMode->NBins;i++){
    TempProf.rstds[i] = Wgt*Profile->rstds[i];
    TempProf.rstdq[i] = Wgt*(cos(AngRot)*Profile->rstdq[i] +
			     sin(AngRot)*Profile->rstdu[i]);
    TempProf.rstdu[i] = Wgt*(cos(AngRot)*Profile->rstdu[i] -
			     sin(AngRot)*Profile->rstdq[i]);
    TempProf.rstdv[i] = Profile->rstdv[i];   
  }
#endif  

  memcpy(Profile,&TempProf,sizeof(struct StdProfs));
  


}




