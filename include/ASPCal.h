/* For .h files and functions particular to ASPCal.c */

// #include "CalCmdLine.h" 

/* Functions that specifically use above structures and CalCmdLine.h */

int GetCalOpt(struct RunVars *, struct CalVars *, 
		 Cmdline *, struct ASPHdr *);
int GetContOpt(struct RunVars *, struct CalVars *, 
		 Cmdline *, struct ASPHdr *);
int GetCalData(struct ASPHdr *, struct SubHdr *, struct RunVars *, 
		  fitsfile *, double **, double **, double **, double **);

int GetPhases(struct ASPHdr *, struct RunVars *,struct CalVars *, 
	      double **, double **, int *, int *);

int GetCalGain(struct CalVars *, 
	       double *, double *, double *, 
	       int, int, 
	       double *, double *, 
	       double *, double *, double *);

int GetCalTsys(struct CalVars *, 
	       double *, double *, double *, 
	       int, int, 
	       double *, double *, 
	       double *, double *, double *);
