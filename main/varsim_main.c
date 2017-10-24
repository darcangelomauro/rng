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

    double init_G, final_G, incr_G;
    int init_dim, final_dim;
    int rep_G, rep_dim;

    printf("Input initial G, final G, and how many values must be simulated in between\n");
    scanf("%lf", &init_G);
    scanf("%lf", &final_G);
    scanf("%d", &rep_G);
    
    // simulation
    multicode_wrapper(P_actionD4D2, P_gamma, incr_G, rep_G, incr_dim, rep_dim, r);
    
    // free random number generator
    gsl_rng_free(r);

}
