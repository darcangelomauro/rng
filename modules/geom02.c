#define GEOM02_C

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



// TYPE (0,2) GEOMETRY

void geom_check02()
{
    if(nH != 0 || nL != 2)
    {
        printf("Error: geometry is not (0,2)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma02()
{

    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gammaL[0], 0, 0, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gammaL[0], 1, 1, gsl_complex_rect(0., -1.));

    // gamma 2
    gsl_matrix_complex_set(gammaL[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaL[1], 1, 0, gsl_complex_rect(-1., 0.));
}

void actionD2t02(double* vecS, int* control)
{
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* L1 = L[1];
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex* L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L1, GSL_COMPLEX_ZERO, L1L1);
    gsl_complex trL1L1 = trace(L1L1);
    gsl_matrix_complex_free(L0L0);
    gsl_matrix_complex_free(L1L1);
    if(control[0])
        vecS[0] = +4*dim*GSL_REAL(trL0L0)+4*dim*GSL_REAL(trL1L1);
    if(control[1])
        vecS[1] = +4*dim*GSL_IMAG(trL0L0)+4*dim*GSL_IMAG(trL1L1);
}

void actionD4D2t02(double* vecS, int* control)
{
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* L1 = L[1];
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex* L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L1, GSL_COMPLEX_ZERO, L1L1);
    gsl_complex trL1L1 = trace(L1L1);
    gsl_matrix_complex* L0L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L1, GSL_COMPLEX_ZERO, L0L1);
    gsl_complex trL0L1 = trace(L0L1);
    gsl_matrix_complex* L1L1L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1L1, L1L1, GSL_COMPLEX_ZERO, L1L1L1L1);
    gsl_complex trL1L1L1L1 = trace(L1L1L1L1);
    gsl_matrix_complex* L0L0L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0L0, L1L1, GSL_COMPLEX_ZERO, L0L0L1L1);
    gsl_complex trL0L0L1L1 = trace(L0L0L1L1);
    gsl_matrix_complex* L0L0L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0L0, L0L0, GSL_COMPLEX_ZERO, L0L0L0L0);
    gsl_complex trL0L0L0L0 = trace(L0L0L0L0);
    gsl_matrix_complex* L0L1L0L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L0L1, L0L1, GSL_COMPLEX_ZERO, L0L1L0L1);
    gsl_complex trL0L1L0L1 = trace(L0L1L0L1);
    gsl_matrix_complex_free(L0L0);
    gsl_matrix_complex_free(L1L1);
    gsl_matrix_complex_free(L0L1);
    gsl_matrix_complex_free(L1L1L1L1);
    gsl_matrix_complex_free(L0L0L1L1);
    gsl_matrix_complex_free(L0L0L0L0);
    gsl_matrix_complex_free(L0L1L0L1);
    if(control[0])
        vecS[0] = G*(+4*dim*GSL_REAL(trL0L0)+4*dim*GSL_REAL(trL1L1)) -8*dim*GSL_REAL(trL0L1L0L1)+16*dim*GSL_REAL(trL0L0L1L1)+16*GSL_REAL(gsl_complex_mul(trL0L1,trL0L1))+12*GSL_REAL(gsl_complex_mul(trL0L0,trL0L0))+4*dim*GSL_REAL(trL1L1L1L1)+4*dim*GSL_REAL(trL0L0L0L0)+8*GSL_REAL(gsl_complex_mul(trL0L0,trL1L1))+12*GSL_REAL(gsl_complex_mul(trL1L1,trL1L1));
    if(control[1])
        vecS[1] = G*(+4*dim*GSL_IMAG(trL0L0)+4*dim*GSL_IMAG(trL1L1)) -8*dim*GSL_IMAG(trL0L1L0L1)+16*dim*GSL_IMAG(trL0L0L1L1)+16*GSL_IMAG(gsl_complex_mul(trL0L1,trL0L1))+12*GSL_IMAG(gsl_complex_mul(trL0L0,trL0L0))+4*dim*GSL_IMAG(trL1L1L1L1)+4*dim*GSL_IMAG(trL0L0L0L0)+8*GSL_IMAG(gsl_complex_mul(trL0L0,trL1L1))+12*GSL_IMAG(gsl_complex_mul(trL1L1,trL1L1));
}




