#define GEOM20_C

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



// TYPE (2,0) GEOMETRY

void geom_check20()
{
    if(nH != 2 || nL != 0)
    {
        printf("Error: geometry is not (2,0)\n");
        exit(EXIT_FAILURE);
    }
}


// build gamma matrices
void init_gamma20()
{

    for(int i=0; i<nHL; i++)
        gammaH[i] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gammaH[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaH[0], 1, 1, gsl_complex_rect(-1., 0.));

    // gamma 2
    gsl_matrix_complex_set(gammaH[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaH[1], 1, 0, gsl_complex_rect(1., 0.));
}
    
void actionD2t20(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* H1 = H[1];
    gsl_complex trH1 = trace(H1);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H1, H1, GSL_COMPLEX_ZERO, H1H1);
    gsl_complex trH1H1 = trace(H1H1);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(H1H1);
    if(control[0])
        vecS[0] = +4*GSL_REAL(gsl_complex_mul(trH1,trH1))+4*dim*GSL_REAL(trH0H0)+4*dim*GSL_REAL(trH1H1)+4*GSL_REAL(gsl_complex_mul(trH0,trH0));
    if(control[1])
        vecS[1] = +4*GSL_IMAG(gsl_complex_mul(trH1,trH1))+4*dim*GSL_IMAG(trH0H0)+4*dim*GSL_IMAG(trH1H1)+4*GSL_IMAG(gsl_complex_mul(trH0,trH0));
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
    if(control[4])
        vecS[4] = GSL_REAL(trH1);
    if(control[5])
        vecS[5] = GSL_REAL(trH1H1);
}

void actionD4D2t20(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* H1 = H[1];
    gsl_complex trH1 = trace(H1);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H1, H1, GSL_COMPLEX_ZERO, H1H1);
    gsl_complex trH1H1 = trace(H1H1);
    gsl_matrix_complex* H0H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H1, GSL_COMPLEX_ZERO, H0H1);
    gsl_complex trH0H1 = trace(H0H1);
    gsl_matrix_complex* H0H0H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0H1, GSL_COMPLEX_ZERO, H0H0H1);
    gsl_complex trH0H0H1 = trace(H0H0H1);
    gsl_matrix_complex* H0H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H1H1, GSL_COMPLEX_ZERO, H0H1H1);
    gsl_complex trH0H1H1 = trace(H0H1H1);
    gsl_matrix_complex* H1H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H1, H1H1, GSL_COMPLEX_ZERO, H1H1H1);
    gsl_complex trH1H1H1 = trace(H1H1H1);
    gsl_matrix_complex* H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0);
    gsl_complex trH0H0H0 = trace(H0H0H0);
    gsl_matrix_complex* H0H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0H0);
    gsl_complex trH0H0H0H0 = trace(H0H0H0H0);
    gsl_matrix_complex* H0H1H0H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0H1, H0H1, GSL_COMPLEX_ZERO, H0H1H0H1);
    gsl_complex trH0H1H0H1 = trace(H0H1H0H1);
    gsl_matrix_complex* H0H0H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, H1H1, GSL_COMPLEX_ZERO, H0H0H1H1);
    gsl_complex trH0H0H1H1 = trace(H0H0H1H1);
    gsl_matrix_complex* H1H1H1H1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H1H1, H1H1, GSL_COMPLEX_ZERO, H1H1H1H1);
    gsl_complex trH1H1H1H1 = trace(H1H1H1H1);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(H1H1);
    gsl_matrix_complex_free(H0H1);
    gsl_matrix_complex_free(H0H0H1);
    gsl_matrix_complex_free(H0H1H1);
    gsl_matrix_complex_free(H1H1H1);
    gsl_matrix_complex_free(H0H0H0);
    gsl_matrix_complex_free(H0H0H0H0);
    gsl_matrix_complex_free(H0H1H0H1);
    gsl_matrix_complex_free(H0H0H1H1);
    gsl_matrix_complex_free(H1H1H1H1);
    if(control[0])
        vecS[0] = G*(+4*GSL_REAL(gsl_complex_mul(trH1,trH1))+4*dim*GSL_REAL(trH0H0)+4*dim*GSL_REAL(trH1H1)+4*GSL_REAL(gsl_complex_mul(trH0,trH0))) +4*dim*GSL_REAL(trH1H1H1H1)+16*GSL_REAL(gsl_complex_mul(trH0H0H1,trH1))+4*dim*GSL_REAL(trH0H0H0H0)+16*dim*GSL_REAL(trH0H0H1H1)+16*GSL_REAL(gsl_complex_mul(trH0,trH0H1H1))+12*GSL_REAL(gsl_complex_mul(trH0H0,trH0H0))+12*GSL_REAL(gsl_complex_mul(trH1H1,trH1H1))+8*GSL_REAL(gsl_complex_mul(trH0H0,trH1H1))+16*GSL_REAL(gsl_complex_mul(trH0H1,trH0H1))+16*GSL_REAL(gsl_complex_mul(trH1,trH1H1H1))-8*dim*GSL_REAL(trH0H1H0H1)+16*GSL_REAL(gsl_complex_mul(trH0,trH0H0H0));
    if(control[1])
        vecS[1] = G*(+4*GSL_IMAG(gsl_complex_mul(trH1,trH1))+4*dim*GSL_IMAG(trH0H0)+4*dim*GSL_IMAG(trH1H1)+4*GSL_IMAG(gsl_complex_mul(trH0,trH0))) +4*dim*GSL_IMAG(trH1H1H1H1)+16*GSL_IMAG(gsl_complex_mul(trH0H0H1,trH1))+4*dim*GSL_IMAG(trH0H0H0H0)+16*dim*GSL_IMAG(trH0H0H1H1)+16*GSL_IMAG(gsl_complex_mul(trH0,trH0H1H1))+12*GSL_IMAG(gsl_complex_mul(trH0H0,trH0H0))+12*GSL_IMAG(gsl_complex_mul(trH1H1,trH1H1))+8*GSL_IMAG(gsl_complex_mul(trH0H0,trH1H1))+16*GSL_IMAG(gsl_complex_mul(trH0H1,trH0H1))+16*GSL_IMAG(gsl_complex_mul(trH1,trH1H1H1))-8*dim*GSL_IMAG(trH0H1H0H1)+16*GSL_IMAG(gsl_complex_mul(trH0,trH0H0H0));
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
    if(control[4])
        vecS[4] = GSL_REAL(trH1);
    if(control[5])
        vecS[5] = GSL_REAL(trH1H1);
}
