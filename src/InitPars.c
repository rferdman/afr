
/* ================================================================ */
/* Fill ASP header with some fake values.. For test purposes.       */
/* INPUT : ASPHdr -- structure of ASP header.                       */
/* R. Ramachandran, 5-November-2003, Berkeley.                        */
/* ================================================================ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ASPHeader.h"

void InitPars(struct ASPHdr *hdr)
{
 /*  int    i; */

  strcpy(hdr->gen.ScanName, "\0");
  strcpy(hdr->gen.SoftVer, "\0");
  strcpy(hdr->gen.Platform, "0");
  strcpy(hdr->gen.CommentOper, "\0");
  strcpy(hdr->gen.HdrVer, "\0");
  strcpy(hdr->gen.FitsType, "\0");
  strcpy(hdr->gen.Observer, "\0");
  strcpy(hdr->gen.ProjID, "\0");
  strcpy(hdr->gen.FEName, "\0");
  strcpy(hdr->gen.FEPol, "\0");
  strcpy(hdr->gen.BEName, "\0");
  strcpy(hdr->gen.BEConfFile, "\0");
  strcpy(hdr->gen.ObsMode, "\0");
  strcpy(hdr->obs.ObsvtyCode, "\0");

  strcpy(hdr->target.PSRName, "\0");
  hdr->target.RA                         = 0.0;
  hdr->target.Dec                        = 0.0;
  hdr->target.Epoch                      = 0.0;
  strcpy(hdr->target.CoordMode, "\0");
  //  hdr->target.StartCrd1                  = 0.0;
  //  hdr->target.StartCrd2                  = 0.0;
  strcpy(hdr->target.StartCrd1, "\0");
  strcpy(hdr->target.StartCrd2, "\0");
  //  hdr->target.StopCrd1                   = 0.0;
  //  hdr->target.StopCrd2                   = 0.0;
  strcpy(hdr->target.StopCrd1, "\0");
  strcpy(hdr->target.StopCrd2, "\0");
  strcpy(hdr->target.TrackMode, "\0");

  hdr->obs.ObsLength                     = 0.0;
  strcpy(hdr->obs.StartDate, "\0");
  strcpy(hdr->obs.StartUT, "\0");
  hdr->obs.IMJDStart                     = 0;
  hdr->obs.StartTime                     = 0;
  hdr->obs.StartFSec                     = 0.0;
  hdr->obs.ClockOffset                   = 0.0;
  hdr->target.StartLST                   = 0.0;
  hdr->obs.IonRM                         = 0.0;
  strcpy(hdr->obs.OPString, "\0");
  hdr->obs.OPScale                       = 0;
  hdr->obs.DynDC                         = 0;
  hdr->obs.DynTime                       = 0.0;
  hdr->obs.NBitPerSamp                   = 0;
  hdr->obs.InterfString                  = 0;
  strcpy(hdr->obs.ZapFileName, "\0");

  hdr->obs.NChan                         = 0;
  hdr->obs.SampInterval                  = 0.0;
  hdr->asp.SplitFac                      = 0;

  hdr->obs.DM                            = 0.0;
  hdr->obs.DMMethod                      = 0;
  hdr->obs.RM                            = 0.0;
  hdr->obs.RMMethod                      = 0;

  hdr->redn.RIMJDFirstSamp              = 0;
  hdr->redn.RfMJDFirstSamp              = 0.0;
  hdr->redn.RNTimeDumps                 = 0;
  hdr->redn.TDump                       = 0.0;
  hdr->redn.RNBinTimeDump               = 0;
  hdr->redn.RChanNum                    = 0;
  hdr->redn.RChanWidth                  = 0.0;
  strcpy(hdr->redn.ROPString, "\0");
  hdr->redn.RDM                         = 0.0;
  hdr->redn.RRM                         = 0.0;
  hdr->redn.RPeriod                     = 0.0;
  hdr->redn.RFoldType                   = 0;
  hdr->redn.RNPolyCoeff                 = 0;
  strcpy(hdr->redn.RDateString, "\0");
  strcpy(hdr->redn.RUserName, "\0");
  hdr->redn.ROPModeCode                 = 0;
  strcpy(hdr->redn.RCommString, "\0");

  return;
}

