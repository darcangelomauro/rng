#define DIRAC_C

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "update.h"
#include "matop.h"
#include "math.h"
#include "global.h"

// builds dirac operator optimizing computations (it's geometry independent)
/*
void build_dirac_fast()
{
    gsl_matrix_complex* a = gsl_matrix_complex_calloc(dimG*dim, dimG*dim);
    gsl_matrix_complex* I = gsl_matrix_complex_alloc(dim, dim);
    gsl_matrix_complex_set_identity(I);
    gsl_matrix_complex_set_zero(DIRAC);

    // first part
    for(int i=0; i<nH; i++)
        tensor_refined(gammaH[i], H[i], a, GSL_COMPLEX_ONE, 1, 1, 1);
    for(int i=0; i<nL; i++)
        tensor_refined(gammaL[i], L[i], a, gsl_complex_rect(0., 1.), 1, 1, 1);

    tensor(a, I, DIRAC);


    // second part
    for(int i=0; i<nH; i++)
        tensor_refined(gammaH1[i], H[i], DIRAC, GSL_COMPLEX_ONE, 1, 0, 1); 
    for(int i=0; i<nL; i++)
        tensor_refined(gammaL1[i], L[i], DIRAC, gsl_complex_rect(0., -1.), 1, 0, 1); 


    gsl_matrix_complex_free(a);
    gsl_matrix_complex_free(I);
}
*/

// builds dirac operator using bruteforce (it's geometry independent)
// (somewhat slower)
void build_dirac()
{
    gsl_matrix_complex* a = gsl_matrix_complex_alloc(dim*dim, dim*dim);
    gsl_matrix_complex* b = gsl_matrix_complex_alloc(dimD, dimD);
    gsl_matrix_complex_set_zero(DIRAC);
    
    for( int i=0; i<nH; i++)
    {
        anticommutator(H[i], a);
        tensor(gammaH[i], a, b);
        gsl_matrix_complex_add(DIRAC, b);
    }
    for( int i=0; i<nL; i++)
    {
        commutator(L[i], a);
        tensor(gammaL[i], a, b);
        gsl_matrix_complex_scale(b, gsl_complex_rect(0., 1.));
        gsl_matrix_complex_add(DIRAC, b);
    }

    gsl_matrix_complex_free(a);
    gsl_matrix_complex_free(b);
}

gsl_complex actionD2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);

    gsl_complex trD2 = trace(D2);

    gsl_matrix_complex_free(D2);

    return trD2;
}

gsl_complex actionD4D2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_matrix_complex* D4 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, D2, D2, GSL_COMPLEX_ZERO, D4);

    gsl_complex trD2 = trace(D2);
    gsl_complex trD4 = trace(D4);

    gsl_matrix_complex_free(D2);
    gsl_matrix_complex_free(D4);

    return gsl_complex_add(gsl_complex_mul_real(trD2, G), trD4);
}


