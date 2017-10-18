#define GEOM01_C

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

// TYPE (0,1) GEOMETRY

void geom_check01()
{
    if(nH != 0 || nL != 1)
    {
        printf("Error: geometry is not (0,1)\n");
        exit(EXIT_FAILURE);
    }
}

void init_gamma01()
{
    gammaL[0] = gsl_matrix_complex_calloc(1, 1);
    gsl_matrix_complex_set(gammaL[0], 0, 0, gsl_complex_rect(0., -1.));
}

void actionD2t01(double* vecS, int* control)
{
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex_free(L0L0);
    if(control[0])
        vecS[0] = +2*dim*GSL_REAL(trL0L0);
    if(control[1])
        vecS[1] = +2*dim*GSL_IMAG(trL0L0);
}

void actionD4D2t01(double* vecS, int* control)
{
    gsl_matrix_complex* L0 = L[0];
    gsl_matrix_complex* L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0, L0, GSL_COMPLEX_ZERO, L0L0);
    gsl_complex trL0L0 = trace(L0L0);
    gsl_matrix_complex* L0L0L0L0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L0L0, L0L0, GSL_COMPLEX_ZERO, L0L0L0L0);
    gsl_complex trL0L0L0L0 = trace(L0L0L0L0);
    gsl_matrix_complex_free(L0L0);
    gsl_matrix_complex_free(L0L0L0L0);
    if(control[0])
        vecS[0] = G*(+2*dim*GSL_REAL(trL0L0)) +6*GSL_REAL(gsl_complex_mul(trL0L0,trL0L0))+2*dim*GSL_REAL(trL0L0L0L0);
    if(control[1])
        vecS[1] = G*(+2*dim*GSL_IMAG(trL0L0)) +6*GSL_IMAG(gsl_complex_mul(trL0L0,trL0L0))+2*dim*GSL_IMAG(trL0L0L0L0);
}
