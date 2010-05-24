/* phase_calc.c 
 * Rewrite of Rob's phase calc.  This one will only 
 * do one pulsar :) Also code is commented now.
 * Uses same Polyco struct.
 * Inputs:  pc = pointer to array of polyco sets for the pulsar.
 *          num_pc = number of polyco sets pointed to by pc.
 *          imjd, fmjd = int, frac part of MJD to calc phase of.
 * Outputs: phase = phase of the pulsar at given MJD. (units 0->1.0)
 *          psrfreq = freq of the pulsar at given MJD. (units???)
 * Returns: Index of polyco set used on success.
 *          -1 on error (MJD out of range of all sets).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "phase_calc.h"

int PhaseCalc(struct Polyco *pc, int num_pc, double imjd, double fmjd, 
        double *phase, double *psrfreq) {

    double dtmin;
    int i, j, pc_set_used=-1;

    *psrfreq=0.0; 
    /* Loop over polyco sets to find applicable one */
    for (j=0; j<num_pc; j++) {
        dtmin = 1440.0*((imjd - pc[j].MjdMidInt)
            + (fmjd - pc[j].MjdMidFrac));   /* Time diff in minutes */
        /* If we're in a valid time range, do the computation */
        if (fabs(dtmin) < pc[j].NMinutes/2.0) {
            *phase = pc[j].Coeff[pc[j].NCoeff-1]; /* Init phase for loop */
            /* Loop over coeffs */
            for (i=pc[j].NCoeff-1; i>0; --i) {
                *psrfreq = dtmin*(*psrfreq) + (double)i*pc[j].Coeff[i];
                *phase = dtmin*(*phase) + pc[j].Coeff[i-1];
            }
            /* Add in DC terms and scale */
            *psrfreq = *psrfreq/60.0 + pc[j].FRotRef;
            *phase += pc[j].PhRotRef + dtmin*60.0*pc[j].FRotRef;
            *phase -= floor(*phase);
            pc_set_used = j;  /* Note which polyco set this was */
            break;
        }
    }

    //printf("")

    /* If no valid polyco set found, return error */
    if (pc_set_used == -1) {
        printf("MJD %9.3f out of valid polyco range.\n", imjd+fmjd);
        return(-1);
    }

    return(pc_set_used);

}

