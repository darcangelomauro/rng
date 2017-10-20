#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "time.h"
#include "global.h"
#include "macros.h"

#define STEP_G 0.2
#define REP_G 20

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    
    // simulation
    for(int i=0; i<REP_G; i++)
    {
        double ar = simulation(P_actionD4D2, 0, P_gamma, r);
        printf("dim: %d,    G: %lf\n", dim, G);
        printf("acceptance rate: %lf\n", ar);
        G -= STEP_G;
        overwrite_init_file();
    }
    
    // free random number generator
    gsl_rng_free(r);

}
