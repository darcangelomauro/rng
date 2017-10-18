#define STATISTICS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "statistics.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include "fileop.h"
#include "matop.h"


/*
void alloc_input_filename(char* data, char* simM, char* code)
{
    int n = strlen(code);
    recover_filename(data, simM, code, n);
}
void alloc_input_filename2(char* evlM, char* evlD, char* asys, char* code)
{
    int n = strlen(code);
    recover_filename2(evlM, evlD, asys, code, n);
}
*/

void read_input(FILE* fsimM, int dim_s, int nH_s, int nL_s, int n, FILE* fout)
{
    gsl_matrix_complex** L_s = malloc(3*sizeof(gsl_matrix_complex*));
    for(int i=0; i<3; i++)
        L_s[i] = gsl_matrix_complex_calloc(dim_s, dim_s);

    double averageRE = 0;
    double averageIM = 0;

    for(int l=0; l<n; l++)
    {
        for(int i=0; i<nH_s; i++)
        {
            for(int j=0; j<dim_s; j++)
            {
                for(int k=0; k<dim_s; k++)
                {
                    double re, im;
                    int read = 0;
                    read += fscanf(fsimM, "%lf", &re);
                    read += fscanf(fsimM, "%lf", &im);
                    if(read != 2)
                    {
                        printf("Error: failed to read from file\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        for(int i=0; i<nL_s; i++)
        {
            for(int j=0; j<dim_s; j++)
            {
                for(int k=0; k<dim_s; k++)
                {
                    double re, im;
                    int read = 0;
                    read += fscanf(fsimM, "%lf", &re);
                    read += fscanf(fsimM, "%lf", &im);
                    if(read != 2)
                    {
                        printf("Error: failed to read from file\n");
                        exit(EXIT_FAILURE);
                    }
                    if(i > 2)
                        gsl_matrix_complex_set(L_s[i-3], j, k, gsl_complex_rect(re, im));
                }
            }
        }
        gsl_matrix_complex* M1 = gsl_matrix_complex_calloc(dim_s, dim_s);
        gsl_matrix_complex* M2 = gsl_matrix_complex_calloc(dim_s, dim_s);
        gsl_matrix_complex* A1 = gsl_matrix_complex_calloc(dim_s, dim_s);
        gsl_matrix_complex* A2 = gsl_matrix_complex_calloc(dim_s, dim_s);
        gsl_matrix_complex* L2 = gsl_matrix_complex_calloc(dim_s, dim_s);
        gsl_matrix_complex* F2 = gsl_matrix_complex_calloc(dim_s, dim_s);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L_s[2], L_s[2], GSL_COMPLEX_ZERO, L2);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L_s[0], L_s[1], GSL_COMPLEX_ZERO, M1);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, L_s[1], L_s[0], GSL_COMPLEX_ZERO, M2);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, M1, L_s[2], GSL_COMPLEX_ZERO, A1);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, M2, L_s[2], GSL_COMPLEX_ZERO, A2);

        gsl_matrix_complex_sub(M1, M2);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, M1, M1, GSL_COMPLEX_ZERO, F2);
        gsl_matrix_complex_sub(A1, A2);

        // at this point:
        // M1 = [L0, L1]
        // A1 = [L0, L1]*L2
        // F2 = M1^2
        // L2 = L2^2
        gsl_complex trA1 = trace(A1);
        gsl_complex res = gsl_complex_mul(trA1, gsl_complex_rect(0., -1.));
        double trF2 = -trace_herm(F2);
        double trL2 = -trace_herm(L2);
        

        // we want to output trA1/(trF2*trL2)

        double DEN = trF2*trL2;
        gsl_complex res_norm = gsl_complex_div_real(res, DEN);

        fprintf(fout, "%.15lf %.15lf\n", GSL_REAL(res_norm), GSL_IMAG(res_norm));

        averageRE += GSL_REAL(res_norm);
        averageIM += GSL_IMAG(res_norm);

        gsl_matrix_complex_free(M1);
        gsl_matrix_complex_free(M2);
        gsl_matrix_complex_free(A1);
        gsl_matrix_complex_free(A2);
        gsl_matrix_complex_free(F2);
        gsl_matrix_complex_free(L2);
    }

    for(int i=0; i<3; i++)
        gsl_matrix_complex_free(L_s[i]);
    free(L_s);

    printf("%lf %lf\n", averageRE/(double)n, averageIM/(double)n);
}

// returns size of binned vector
int binned_size(int n, int dimbin)
{
    int size = n/dimbin;
    if((double)n/dimbin != size)
        size++;

    return size;
}


// bins a n components vector with bins of size dimbin
// and returns the binned vector
double* binned_vector(double* v, int n, int dimbin)
{
    if(dimbin > n)
    {
        printf("Error, bin size cannot be greater than vector length\n");
        exit(EXIT_FAILURE);
    }

    int size = binned_size(n, dimbin);
    double* u = calloc(size, sizeof(double));
    int* count = calloc(size, sizeof(double));

    for(int i=0; i<n; i++)
    {
        int idx = i/dimbin;
        if(idx >= size)
        {
            printf("Error index out of range\n");
            exit(EXIT_FAILURE);
        }

        u[idx] += v[i];
        count[idx]++;
    }

    for(int i=0; i<size; i++)
        u[i] /= (double)count[i];

    return u;
}






