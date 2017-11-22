#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include "time.h"
#include "global.h"
#include "macros.h"
#include "matop.h"

#define REP 500
#define ACTION_B actionD4D2_bruteforce
#define ACTION actionD4D2t22

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // file
    FILE* fS = fopen("D4D2t22.txt", "w");
    FILE* f2 = fopen("D4D2t22_time.txt", "w");

    // action vector
    double vecS[14];
    int control[14] = {1,1,0,0,0,0,0,0,0,0,0,0,0,0};
    
    init_data();
    GEOM_CHECK();
    init_hot(ACTION, P_gamma, r);
    
    clock_t start1;
    clock_t cumul1 = 0;
    clock_t start2;
    clock_t cumul2 = 0;
    // simulation
    for(int i=0; i<REP; i++)
    {
        for(int k=0; k<nH; k++)
            generate_HL(H[k], 0, dim, r);  


        vecS[0] = 0.;
        vecS[1] = 0.;
        
        start1 = clock();
        ACTION(vecS, control);
        cumul1 += clock()-start1;

        start2 = clock();
        gsl_complex cS = ACTION_B();
        cumul2 += clock()-start2;
        
        fprintf(fS, "%.15lf %.15lf %.15lf %.15lf\n", vecS[0], GSL_REAL(cS), vecS[1], GSL_IMAG(cS));

    }

    fprintf(f2, "nH = %d, nL = %d, REP = %d, n = %d:\n", nH, nL, REP, dim);
    fprintf(f2, "decomp time %lf\n", (double)cumul1/(double)CLOCKS_PER_SEC);
    fprintf(f2, "brutef time %lf\n", (double)cumul2/(double)CLOCKS_PER_SEC);

    // free random number generator
    gsl_rng_free(r);
    fclose(fS);
    fclose(f2);
}
