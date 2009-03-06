/* 

   Reads input data file header info corresponding to each dump and 
   calcuates output dump/channel info.
   
   Specifically, the information is:
   
   DumpMiddleSecs
   DumpRefPhase[channel]
   DumpRefPerdios[channel]
   
   from the SubHDr structure.  -- RDF, 09/09/2008 
   
*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ASPCommon.h"
#include "ASPFitsReader.h"

int MakeSubHDr( struct RunVars *RunMode, , Cmdline *Cmd, struct ASPHdr *hdr, struct SubHdr *SubInHdr, struct SubHdr *SubOutHdr )
{
  
  /* Assuming FITS file is open.  We already know NDumps.  */

  










  return 0;

}
