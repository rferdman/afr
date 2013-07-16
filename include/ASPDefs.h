#define STRINGLEN  32
//#define NCHMAX   128   /* Max no. of DFB channels */ 
//#define NCHMAX   2048    /* Max no. of DFB channels */ 
#define NDSMAX    4    /* Max no. of data servers */ 
#ifndef Boolean        /* Define Boolean type definition */ 
#define Boolean int   
#endif 
#ifndef FALSE          /* also defined in "utils.h" */ 
#define FALSE   0 
#endif 
#ifndef TRUE           /* also defined in "utils.h" */ 
#define TRUE    1 
#endif 

#define NBINMAX 4096
#define NB      128
#define MAXSPLIT 32
#define NFILEMAX 1000
#define MAXDUMPS 1024
#define MAXOMIT  524288
#define NCHMAX 2048
#define NCOLMAX 999

//#define TWOPI 6.2831853071795864
/* MORE DECIMALS...!  */
#define TWOPI 6.2831853071795864769252867665590057683943387987502
#define PI    3.1415926535897932384626433832795028841971693993751

#define DFFAC 2.410e-4

/* THRESHOLD WHEN COMPARING TWO DOUBLES */
//#define DBLEPS 0.000001
#define DBLEPS 0.00011
/* Tolerance for determining zero values */
#define FLTTOL 0.0000001


