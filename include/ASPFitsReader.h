#include "CmdLine.h" 

int  GetOptions(struct RunVars *, struct CalVars *, 
		Cmdline *, struct ASPHdr *);
int  GetOmit(struct ASPHdr *, Cmdline *, struct RunVars *);
int  GetChans(struct ASPHdr *, Cmdline *, struct RunVars *);
int  ReadCal(struct ASPHdr *, struct RunVars *, struct CalVars *, 
	    double **);
void PrintLog(struct RunVars *, struct ASPHdr *, Cmdline *);
// int MakePoly(Cmdline *, struct ASPHdr *);
