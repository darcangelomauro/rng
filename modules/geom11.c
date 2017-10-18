#define GEOM11_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include "geom.h"
#include "matop.h"
#include "global.h"



// TYPE (1,1) GEOMETRY

void geom_check11()
{
    if(nH != 1 || nL != 1)
    {
        printf("Error: geometry is not (1,1)\n");
        exit(EXIT_FAILURE);
    }
}


// build gamma matrices
void init_gamma11()
{

    gammaH[0] = gsl_matrix_complex_calloc(2, 2);
    gammaL[0] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gammaH[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaH[0], 1, 1, gsl_complex_rect(-1., 0.));

    // gamma 2
    gsl_matrix_complex_set(gammaL[0], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaL[0], 1, 0, gsl_complex_rect(-1., 0.));
}

void actionD2t11(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(L0L0);
    if(control[0])
        vecS[0] = +4*GSL_REAL(gsl_complex_mul(trH0,trH0))+4*dim*GSL_REAL(trL0L0)+4*dim*GSL_REAL(trH0H0);
    if(control[1])
        vecS[1] = +4*GSL_IMAG(gsl_complex_mul(trH0,trH0))+4*dim*GSL_IMAG(trL0L0)+4*dim*GSL_IMAG(trH0H0);
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}

void actionD4D2t11(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* H0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L0, GSL_COMPLEX_ZERO, H0L0);
    gsl_complex trH0L0 = trace(H0L0);
    gsl_matrix_complex* H0L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L0L0, GSL_COMPLEX_ZERO, H0L0L0);
    gsl_complex trH0L0L0 = trace(H0L0L0);
    gsl_matrix_complex* H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0);
    gsl_complex trH0H0H0 = trace(H0H0H0);
    gsl_matrix_complex* H0H0L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, L0L0, GSL_COMPLEX_ZERO, H0H0L0L0);
    gsl_complex trH0H0L0L0 = trace(H0H0L0L0);
    gsl_matrix_complex* H0L0H0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L0, H0L0, GSL_COMPLEX_ZERO, H0L0H0L0);
    gsl_complex trH0L0H0L0 = trace(H0L0H0L0);
    gsl_matrix_complex* L0L0L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0L0, L0L0, GSL_COMPLEX_ZERO, L0L0L0L0);
    gsl_complex trL0L0L0L0 = trace(L0L0L0L0);
    gsl_matrix_complex* H0H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0H0);
    gsl_complex trH0H0H0H0 = trace(H0H0H0H0);
    gsl_matrix_complex_free(L0L0);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(H0L0);
    gsl_matrix_complex_free(H0L0L0);
    gsl_matrix_complex_free(H0H0H0);
    gsl_matrix_complex_free(H0H0L0L0);
    gsl_matrix_complex_free(H0L0H0L0);
    gsl_matrix_complex_free(L0L0L0L0);
    gsl_matrix_complex_free(H0H0H0H0);
    if(control[0])
        vecS[0] = G*(+4*dim*GSL_REAL(trH0H0)+4*GSL_REAL(gsl_complex_mul(trH0,trH0))+4*dim*GSL_REAL(trL0L0)) -8*dim*GSL_REAL(trH0L0H0L0)+16*dim*GSL_REAL(trH0H0L0L0)+4*dim*GSL_REAL(trL0L0L0L0)-16*GSL_REAL(gsl_complex_mul(trH0L0,trH0L0))+12*GSL_REAL(gsl_complex_mul(trH0H0,trH0H0))+4*dim*GSL_REAL(trH0H0H0H0)+16*GSL_REAL(gsl_complex_mul(trH0,trH0L0L0))+16*GSL_REAL(gsl_complex_mul(trH0,trH0H0H0))+8*GSL_REAL(gsl_complex_mul(trH0H0,trL0L0))+12*GSL_REAL(gsl_complex_mul(trL0L0,trL0L0));
    if(control[1])
        vecS[1] = G*(+4*dim*GSL_IMAG(trH0H0)+4*GSL_IMAG(gsl_complex_mul(trH0,trH0))+4*dim*GSL_IMAG(trL0L0)) -8*dim*GSL_IMAG(trH0L0H0L0)+16*dim*GSL_IMAG(trH0H0L0L0)+4*dim*GSL_IMAG(trL0L0L0L0)-16*GSL_IMAG(gsl_complex_mul(trH0L0,trH0L0))+12*GSL_IMAG(gsl_complex_mul(trH0H0,trH0H0))+4*dim*GSL_IMAG(trH0H0H0H0)+16*GSL_IMAG(gsl_complex_mul(trH0,trH0L0L0))+16*GSL_IMAG(gsl_complex_mul(trH0,trH0H0H0))+8*GSL_IMAG(gsl_complex_mul(trH0H0,trL0L0))+12*GSL_IMAG(gsl_complex_mul(trL0L0,trL0L0));
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}
