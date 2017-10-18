#define GEOM03_C

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

// TYPE (0,3) GEOMETRY

void geom_check03()
{
    if(nH != 1 || nL != 3)
    {
        printf("Error: geometry is not (0,3)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma03()
{

    gammaH[0] = gsl_matrix_complex_calloc(2, 2);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(2, 2);

    // gamma 0
    gsl_matrix_complex_set(gammaH[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaH[0], 1, 1, gsl_complex_rect(1., 0.));

    // gamma 1
    gsl_matrix_complex_set(gammaL[0], 0, 0, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gammaL[0], 1, 1, gsl_complex_rect(0., -1.));

    // gamma 2
    gsl_matrix_complex_set(gammaL[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gammaL[1], 1, 0, gsl_complex_rect(-1., 0.));
    
    // gamma 3
    gsl_matrix_complex_set(gammaL[2], 0, 1, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gammaL[2], 1, 0, gsl_complex_rect(0., 1.));
}

void actionD2t03(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* L1 = L[0];
    gsl_matrix_complex* L2 = L[1];
    gsl_matrix_complex* L3 = L[2];
    gsl_matrix_complex* L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L1, GSL_COMPLEX_ZERO, L1L1);
    gsl_complex trL1L1 = trace(L1L1);
    gsl_matrix_complex* L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2, L2, GSL_COMPLEX_ZERO, L2L2);
    gsl_complex trL2L2 = trace(L2L2);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L3, L3, GSL_COMPLEX_ZERO, L3L3);
    gsl_complex trL3L3 = trace(L3L3);
    gsl_matrix_complex_free(L1L1);
    gsl_matrix_complex_free(L2L2);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(L3L3);
    if(control[0])
        vecS[0] = +4*dim*GSL_REAL(trL3L3)+4*dim*GSL_REAL(trL2L2)+4*dim*GSL_REAL(trH0H0)+4*dim*GSL_REAL(trL1L1)+4*GSL_REAL(gsl_complex_mul(trH0,trH0));
    if(control[1])
        vecS[1] = +4*dim*GSL_IMAG(trL3L3)+4*dim*GSL_IMAG(trL2L2)+4*dim*GSL_IMAG(trH0H0)+4*dim*GSL_IMAG(trL1L1)+4*GSL_IMAG(gsl_complex_mul(trH0,trH0));
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}

void actionD4D2t03(double* vecS, int* control)
{
    gsl_matrix_complex* H0 = H[0];
    gsl_complex trH0 = trace(H0);
    gsl_matrix_complex* L1 = L[0];
    gsl_matrix_complex* L2 = L[1];
    gsl_matrix_complex* L3 = L[2];
    gsl_matrix_complex* H0L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L1, GSL_COMPLEX_ZERO, H0L1);
    gsl_complex trH0L1 = trace(H0L1);
    gsl_matrix_complex* H0L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L3, GSL_COMPLEX_ZERO, H0L3);
    gsl_complex trH0L3 = trace(H0L3);
    gsl_matrix_complex* L1L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L2, GSL_COMPLEX_ZERO, L1L2);
    gsl_complex trL1L2 = trace(L1L2);
    gsl_matrix_complex* H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0, GSL_COMPLEX_ZERO, H0H0);
    gsl_complex trH0H0 = trace(H0H0);
    gsl_matrix_complex* L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2, L2, GSL_COMPLEX_ZERO, L2L2);
    gsl_complex trL2L2 = trace(L2L2);
    gsl_matrix_complex* H0L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L2, GSL_COMPLEX_ZERO, H0L2);
    gsl_complex trH0L2 = trace(H0L2);
    gsl_matrix_complex* L2L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2, L3, GSL_COMPLEX_ZERO, L2L3);
    gsl_complex trL2L3 = trace(L2L3);
    gsl_matrix_complex* L1L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L3, GSL_COMPLEX_ZERO, L1L3);
    gsl_complex trL1L3 = trace(L1L3);
    gsl_matrix_complex* L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L1, GSL_COMPLEX_ZERO, L1L1);
    gsl_complex trL1L1 = trace(L1L1);
    gsl_matrix_complex* L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L3, L3, GSL_COMPLEX_ZERO, L3L3);
    gsl_complex trL3L3 = trace(L3L3);
    gsl_matrix_complex* L3L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L3, L2, GSL_COMPLEX_ZERO, L3L2);
    gsl_matrix_complex* L3L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L3, L1, GSL_COMPLEX_ZERO, L3L1);
    gsl_matrix_complex* L2L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2, L1, GSL_COMPLEX_ZERO, L2L1);
    gsl_matrix_complex* H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0);
    gsl_complex trH0H0H0 = trace(H0H0H0);
    gsl_matrix_complex* L1L3L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L3L2, GSL_COMPLEX_ZERO, L1L3L2);
    gsl_complex trL1L3L2 = trace(L1L3L2);
    gsl_matrix_complex* H0L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L1L1, GSL_COMPLEX_ZERO, H0L1L1);
    gsl_complex trH0L1L1 = trace(H0L1L1);
    gsl_matrix_complex* H0L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L3L3, GSL_COMPLEX_ZERO, H0L3L3);
    gsl_complex trH0L3L3 = trace(H0L3L3);
    gsl_matrix_complex* L1L2L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1, L2L3, GSL_COMPLEX_ZERO, L1L2L3);
    gsl_complex trL1L2L3 = trace(L1L2L3);
    gsl_matrix_complex* H0L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0, L2L2, GSL_COMPLEX_ZERO, H0L2L2);
    gsl_complex trH0L2L2 = trace(H0L2L2);
    gsl_matrix_complex* H0L2L3L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L2, L3L1, GSL_COMPLEX_ZERO, H0L2L3L1);
    gsl_complex trH0L2L3L1 = trace(H0L2L3L1);
    gsl_matrix_complex* H0L3L1L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L3, L1L2, GSL_COMPLEX_ZERO, H0L3L1L2);
    gsl_complex trH0L3L1L2 = trace(H0L3L1L2);
    gsl_matrix_complex* L2L2L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2L2, L3L3, GSL_COMPLEX_ZERO, L2L2L3L3);
    gsl_complex trL2L2L3L3 = trace(L2L2L3L3);
    gsl_matrix_complex* L3L3L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L3L3, L3L3, GSL_COMPLEX_ZERO, L3L3L3L3);
    gsl_complex trL3L3L3L3 = trace(L3L3L3L3);
    gsl_matrix_complex* H0H0L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, L1L1, GSL_COMPLEX_ZERO, H0H0L1L1);
    gsl_complex trH0H0L1L1 = trace(H0H0L1L1);
    gsl_matrix_complex* H0H0L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, L3L3, GSL_COMPLEX_ZERO, H0H0L3L3);
    gsl_complex trH0H0L3L3 = trace(H0H0L3L3);
    gsl_matrix_complex* L1L1L3L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1L1, L3L3, GSL_COMPLEX_ZERO, L1L1L3L3);
    gsl_complex trL1L1L3L3 = trace(L1L1L3L3);
    gsl_matrix_complex* L1L3L1L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L1L3, L1L3, GSL_COMPLEX_ZERO, L1L3L1L3);
    gsl_complex trL1L3L1L3 = trace(L1L3L1L3);
    gsl_matrix_complex* L2L2L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L2L2, L2L2, GSL_COMPLEX_ZERO, L2L2L2L2);
    gsl_complex trL2L2L2L2 = trace(L2L2L2L2);
    gsl_matrix_complex* H0L2L1L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L2, L1L3, GSL_COMPLEX_ZERO, H0L2L1L3);
    gsl_complex trH0L2L1L3 = trace(H0L2L1L3);
    gsl_matrix_complex* H0H0H0H0 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, H0H0, GSL_COMPLEX_ZERO, H0H0H0H0);
    gsl_complex trH0H0H0H0 = trace(H0H0H0H0);
    gsl_matrix_complex* H0L1L3L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L1, L3L2, GSL_COMPLEX_ZERO, H0L1L3L2);
    gsl_complex trH0L1L3L2 = trace(H0L1L3L2);
    gsl_matrix_complex* L1L1L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1L1, L2L2, GSL_COMPLEX_ZERO, L1L1L2L2);
    gsl_complex trL1L1L2L2 = trace(L1L1L2L2);
    gsl_matrix_complex* L1L1L1L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, L1L1, L1L1, GSL_COMPLEX_ZERO, L1L1L1L1);
    gsl_complex trL1L1L1L1 = trace(L1L1L1L1);
    gsl_matrix_complex* H0L1H0L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L1, H0L1, GSL_COMPLEX_ZERO, H0L1H0L1);
    gsl_complex trH0L1H0L1 = trace(H0L1H0L1);
    gsl_matrix_complex* H0H0L2L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H0H0, L2L2, GSL_COMPLEX_ZERO, H0H0L2L2);
    gsl_complex trH0H0L2L2 = trace(H0H0L2L2);
    gsl_matrix_complex* H0L1L2L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L1, L2L3, GSL_COMPLEX_ZERO, H0L1L2L3);
    gsl_complex trH0L1L2L3 = trace(H0L1L2L3);
    gsl_matrix_complex* H0L2H0L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L2, H0L2, GSL_COMPLEX_ZERO, H0L2H0L2);
    gsl_complex trH0L2H0L2 = trace(H0L2H0L2);
    gsl_matrix_complex* H0L3L2L1 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L3, L2L1, GSL_COMPLEX_ZERO, H0L3L2L1);
    gsl_complex trH0L3L2L1 = trace(H0L3L2L1);
    gsl_matrix_complex* L2L3L2L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L2L3, L2L3, GSL_COMPLEX_ZERO, L2L3L2L3);
    gsl_complex trL2L3L2L3 = trace(L2L3L2L3);
    gsl_matrix_complex* L1L2L1L2 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L1L2, L1L2, GSL_COMPLEX_ZERO, L1L2L1L2);
    gsl_complex trL1L2L1L2 = trace(L1L2L1L2);
    gsl_matrix_complex* H0L3H0L3 = gsl_matrix_complex_calloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H0L3, H0L3, GSL_COMPLEX_ZERO, H0L3H0L3);
    gsl_complex trH0L3H0L3 = trace(H0L3H0L3);
    gsl_matrix_complex_free(H0L1);
    gsl_matrix_complex_free(H0L3);
    gsl_matrix_complex_free(L1L2);
    gsl_matrix_complex_free(H0H0);
    gsl_matrix_complex_free(L2L2);
    gsl_matrix_complex_free(H0L2);
    gsl_matrix_complex_free(L2L3);
    gsl_matrix_complex_free(L1L3);
    gsl_matrix_complex_free(L1L1);
    gsl_matrix_complex_free(L3L3);
    gsl_matrix_complex_free(L3L2);
    gsl_matrix_complex_free(L3L1);
    gsl_matrix_complex_free(L2L1);
    gsl_matrix_complex_free(H0H0H0);
    gsl_matrix_complex_free(L1L3L2);
    gsl_matrix_complex_free(H0L1L1);
    gsl_matrix_complex_free(H0L3L3);
    gsl_matrix_complex_free(L1L2L3);
    gsl_matrix_complex_free(H0L2L2);
    gsl_matrix_complex_free(H0L2L3L1);
    gsl_matrix_complex_free(H0L3L1L2);
    gsl_matrix_complex_free(L2L2L3L3);
    gsl_matrix_complex_free(L3L3L3L3);
    gsl_matrix_complex_free(H0H0L1L1);
    gsl_matrix_complex_free(H0H0L3L3);
    gsl_matrix_complex_free(L1L1L3L3);
    gsl_matrix_complex_free(L1L3L1L3);
    gsl_matrix_complex_free(L2L2L2L2);
    gsl_matrix_complex_free(H0L2L1L3);
    gsl_matrix_complex_free(H0H0H0H0);
    gsl_matrix_complex_free(H0L1L3L2);
    gsl_matrix_complex_free(L1L1L2L2);
    gsl_matrix_complex_free(L1L1L1L1);
    gsl_matrix_complex_free(H0L1H0L1);
    gsl_matrix_complex_free(H0H0L2L2);
    gsl_matrix_complex_free(H0L1L2L3);
    gsl_matrix_complex_free(H0L2H0L2);
    gsl_matrix_complex_free(H0L3L2L1);
    gsl_matrix_complex_free(L2L3L2L3);
    gsl_matrix_complex_free(L1L2L1L2);
    gsl_matrix_complex_free(H0L3H0L3);
    if(control[0])
        vecS[0] = G*(+4*dim*GSL_REAL(trL3L3)+4*dim*GSL_REAL(trL2L2)+4*dim*GSL_REAL(trH0H0)+4*dim*GSL_REAL(trL1L1)+4*GSL_REAL(gsl_complex_mul(trH0,trH0))) +8*GSL_REAL(gsl_complex_mul(trL1L1,trL3L3))-16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L2L3L1))+16*dim*GSL_REAL(trH0H0L2L2)-16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L3L1L2))+16*dim*GSL_REAL(trL2L2L3L3)+4*dim*GSL_REAL(trL1L1L1L1)+16*dim*GSL_REAL(trH0H0L1L1)+4*dim*GSL_REAL(trH0H0H0H0)+16*GSL_REAL(gsl_complex_mul(trL1L2,trL1L2))+8*dim*GSL_REAL(trH0L1H0L1)+24*GSL_REAL(gsl_complex_mul(trH0H0,trL2L2))+16*dim*GSL_REAL(trL1L1L2L2)-8*dim*GSL_REAL(trL1L2L1L2)+48*GSL_REAL(gsl_complex_mul(trH0,trH0L3L3))+24*GSL_REAL(gsl_complex_mul(trH0H0,trL1L1))+16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L2L1L3))+48*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(trH0,trL1L3L2)))+24*GSL_REAL(gsl_complex_mul(trH0H0,trL3L3))+16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L3L2L1))-48*GSL_REAL(gsl_complex_mul(trH0L2,trH0L2))+8*dim*GSL_REAL(trH0L3H0L3)+48*GSL_REAL(gsl_complex_mul(trH0,trH0L1L1))-48*GSL_REAL(gsl_complex_mul(trH0L3,trH0L3))+48*GSL_REAL(gsl_complex_mul(trH0,trH0L2L2))+16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L1L3L2))+12*GSL_REAL(gsl_complex_mul(trL3L3,trL3L3))-8*dim*GSL_REAL(trL1L3L1L3)+16*dim*GSL_REAL(trH0H0L3L3)+12*GSL_REAL(gsl_complex_mul(trL1L1,trL1L1))+4*dim*GSL_REAL(trL2L2L2L2)+8*dim*GSL_REAL(trH0L2H0L2)+16*GSL_REAL(gsl_complex_mul(trL2L3,trL2L3))-16*dim*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L1L2L3))-8*dim*GSL_REAL(trL2L3L2L3)+12*GSL_REAL(gsl_complex_mul(trH0H0,trH0H0))+16*GSL_REAL(gsl_complex_mul(trH0,trH0H0H0))+16*dim*GSL_REAL(trL1L1L3L3)-48*GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(trH0,trL1L2L3)))+16*GSL_REAL(gsl_complex_mul(trL1L3,trL1L3))+4*dim*GSL_REAL(trL3L3L3L3)+12*GSL_REAL(gsl_complex_mul(trL2L2,trL2L2))-48*GSL_REAL(gsl_complex_mul(trH0L1,trH0L1))+8*GSL_REAL(gsl_complex_mul(trL2L2,trL3L3))+8*GSL_REAL(gsl_complex_mul(trL1L1,trL2L2));
    if(control[1])
        vecS[1] = G*(+4*dim*GSL_IMAG(trL3L3)+4*dim*GSL_IMAG(trL2L2)+4*dim*GSL_IMAG(trH0H0)+4*dim*GSL_IMAG(trL1L1)+4*GSL_IMAG(gsl_complex_mul(trH0,trH0))) +8*GSL_IMAG(gsl_complex_mul(trL1L1,trL3L3))-16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L2L3L1))+16*dim*GSL_IMAG(trH0H0L2L2)-16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L3L1L2))+16*dim*GSL_IMAG(trL2L2L3L3)+4*dim*GSL_IMAG(trL1L1L1L1)+16*dim*GSL_IMAG(trH0H0L1L1)+4*dim*GSL_IMAG(trH0H0H0H0)+16*GSL_IMAG(gsl_complex_mul(trL1L2,trL1L2))+8*dim*GSL_IMAG(trH0L1H0L1)+24*GSL_IMAG(gsl_complex_mul(trH0H0,trL2L2))+16*dim*GSL_IMAG(trL1L1L2L2)-8*dim*GSL_IMAG(trL1L2L1L2)+48*GSL_IMAG(gsl_complex_mul(trH0,trH0L3L3))+24*GSL_IMAG(gsl_complex_mul(trH0H0,trL1L1))+16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L2L1L3))+48*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(trH0,trL1L3L2)))+24*GSL_IMAG(gsl_complex_mul(trH0H0,trL3L3))+16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L3L2L1))-48*GSL_IMAG(gsl_complex_mul(trH0L2,trH0L2))+8*dim*GSL_IMAG(trH0L3H0L3)+48*GSL_IMAG(gsl_complex_mul(trH0,trH0L1L1))-48*GSL_IMAG(gsl_complex_mul(trH0L3,trH0L3))+48*GSL_IMAG(gsl_complex_mul(trH0,trH0L2L2))+16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L1L3L2))+12*GSL_IMAG(gsl_complex_mul(trL3L3,trL3L3))-8*dim*GSL_IMAG(trL1L3L1L3)+16*dim*GSL_IMAG(trH0H0L3L3)+12*GSL_IMAG(gsl_complex_mul(trL1L1,trL1L1))+4*dim*GSL_IMAG(trL2L2L2L2)+8*dim*GSL_IMAG(trH0L2H0L2)+16*GSL_IMAG(gsl_complex_mul(trL2L3,trL2L3))-16*dim*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), trH0L1L2L3))-8*dim*GSL_IMAG(trL2L3L2L3)+12*GSL_IMAG(gsl_complex_mul(trH0H0,trH0H0))+16*GSL_IMAG(gsl_complex_mul(trH0,trH0H0H0))+16*dim*GSL_IMAG(trL1L1L3L3)-48*GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(trH0,trL1L2L3)))+16*GSL_IMAG(gsl_complex_mul(trL1L3,trL1L3))+4*dim*GSL_IMAG(trL3L3L3L3)+12*GSL_IMAG(gsl_complex_mul(trL2L2,trL2L2))-48*GSL_IMAG(gsl_complex_mul(trH0L1,trH0L1))+8*GSL_IMAG(gsl_complex_mul(trL2L2,trL3L3))+8*GSL_IMAG(gsl_complex_mul(trL1L1,trL2L2));
    if(control[2])
        vecS[2] = GSL_REAL(trH0);
    if(control[3])
        vecS[3] = GSL_REAL(trH0H0);
}

