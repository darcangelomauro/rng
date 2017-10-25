#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "time.h"
#include "global.h"
#include "macros.h"

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    
    // simulation
    simulation(P_actionD2, 0, P_gamma, r);
    
    // free random number generator
    gsl_rng_free(r);
}
