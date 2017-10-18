#define GEOM10_C

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

// TYPE (1,0) GEOMETRY

void geom_check10()
{
    if(nH != 1 || nL != 0)
    {
        printf("Error: geometry is not (1,0)\n");
        exit(EXIT_FAILURE);
    }
}

void init_gamma10()
{
    gammaH[0] = gsl_matrix_complex_calloc(1, 1);
    gsl_matrix_complex_set(gammaH[0], 0, 0, gsl_complex_rect(1., 0.));
}

void actionD2t10(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex_free(H0H0);
    if(control[0])
        vecS[0] = +2*GSL_REAL(gsl_complex_mul(trH0,trH0))+2*dim*GSL_REAL(trH0H0);
    if(control[1])
        vecS[1] = +2*GSL_IMAG(gsl_complex_mul(trH0,trH0))+2*dim*GSL_IMAG(trH0H0);
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}

void actionD4D2t10(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0);
    gsl_complex trH0H0H0 = trace(H0H0H0);
    gsl_matrix_complex* H0H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0H0);
    gsl_complex trH0H0H0H0 = trace(H0H0H0H0);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(H0H0H0);
    gsl_matrix_complex_free(H0H0H0H0);
    if(control[0])
        vecS[0] = G*(+2*GSL_REAL(gsl_complex_mul(trH0,trH0))+2*dim*GSL_REAL(trH0H0)) +8*GSL_REAL(gsl_complex_mul(trH0,trH0H0H0))+6*GSL_REAL(gsl_complex_mul(trH0H0,trH0H0))+2*dim*GSL_REAL(trH0H0H0H0);
    if(control[1])
        vecS[1] = G*(+2*GSL_IMAG(gsl_complex_mul(trH0,trH0))+2*dim*GSL_IMAG(trH0H0)) +8*GSL_IMAG(gsl_complex_mul(trH0,trH0H0H0))+6*GSL_IMAG(gsl_complex_mul(trH0H0,trH0H0))+2*dim*GSL_IMAG(trH0H0H0H0);
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}
