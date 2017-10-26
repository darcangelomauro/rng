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

    double incr_G;
    int incr_dim, rep_G, rep_dim;
    int renorm;
    int r1 = 0;

    printf("The program allows to simulate automatically a range of values of G and dim.\n");
    printf("The outemost loop increases dim, the innermost increases G.\n");
    printf("Input increase in G at each step (can be negative, resulting in a decrease)\n");
    r1 += scanf("%lf", &incr_G);
    printf("Specify how many values of G must be simulated\n");
    r1 += scanf("%d", &rep_G);
    printf("Specify increase in dim at each step (can be negative, resulting in a decrease, but must be integer)\n");
    r1 += scanf("%d", &incr_dim);
    printf("Specify how many values of dim must be simulated\n");
    r1 += scanf("%d", &rep_dim);
    printf("Specify if renormalization on (1) or off (0) \n");
    r1 += scanf("%d", &renorm);
    
    if(r1 != 5)
    {
        printf("Error: failed to read parameters\n");
        exit(EXIT_FAILURE);
    }
    
    // simulation
    multicode_wrapper(P_actionD4D2, P_gamma, renorm, incr_G, rep_G, incr_dim, rep_dim, r);
    
    // free random number generator
    gsl_rng_free(r);
}
