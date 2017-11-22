#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "matop.h"


int main()
{
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(6, 6);
    gsl_matrix_complex* B = gsl_matrix_complex_alloc(3, 3);

    gsl_matrix_complex_set(A, 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(A, 0, 1, gsl_complex_rect(23423423432., 234324324.));
    gsl_matrix_complex_set(A, 0, 2, gsl_complex_rect(3., 2.));
    gsl_matrix_complex_set(A, 0, 3, gsl_complex_rect(4., 3253525.));
    gsl_matrix_complex_set(A, 0, 4, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 0, 5, gsl_complex_rect(4., 5345435.));
    gsl_matrix_complex_set(A, 1, 0, gsl_complex_rect(23432432432., 234324324.));
    gsl_matrix_complex_set(A, 1, 1, gsl_complex_rect(2., 0.));
    gsl_matrix_complex_set(A, 1, 2, gsl_complex_rect(1., 143543543.));
    gsl_matrix_complex_set(A, 1, 3, gsl_complex_rect(2., 2.));
    gsl_matrix_complex_set(A, 1, 4, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 1, 5, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 2, 0, gsl_complex_rect(3., 3.));
    gsl_matrix_complex_set(A, 2, 1, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 2, 2, gsl_complex_rect(5., 0.));
    gsl_matrix_complex_set(A, 2, 3, gsl_complex_rect(23432432432., 23423432432.));
    gsl_matrix_complex_set(A, 2, 4, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 2, 5, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 3, 0, gsl_complex_rect(1., 1.));
    gsl_matrix_complex_set(A, 3, 1, gsl_complex_rect(2., 2.));
    gsl_matrix_complex_set(A, 3, 2, gsl_complex_rect(23423424324., 32432424.));
    gsl_matrix_complex_set(A, 3, 3, gsl_complex_rect(6., 0.));
    gsl_matrix_complex_set(A, 3, 4, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 3, 5, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(A, 4, 0, gsl_complex_rect(1., 1.));
    gsl_matrix_complex_set(A, 4, 1, gsl_complex_rect(2., 2.));
    gsl_matrix_complex_set(A, 4, 2, gsl_complex_rect(2., 3.));
    gsl_matrix_complex_set(A, 4, 3, gsl_complex_rect(6., 0.));
    gsl_matrix_complex_set(A, 4, 4, gsl_complex_rect(4., 0.));
    gsl_matrix_complex_set(A, 4, 5, gsl_complex_rect(434534534., 4435345345.));
    gsl_matrix_complex_set(A, 5, 0, gsl_complex_rect(1., 1.));
    gsl_matrix_complex_set(A, 5, 1, gsl_complex_rect(2., 2.));
    gsl_matrix_complex_set(A, 5, 2, gsl_complex_rect(2., 3.));
    gsl_matrix_complex_set(A, 5, 3, gsl_complex_rect(6., 0.));
    gsl_matrix_complex_set(A, 5, 4, gsl_complex_rect(435345345., 345345345.));
    gsl_matrix_complex_set(A, 5, 5, gsl_complex_rect(4., 0.));
    printmat_complex(A, 6, 6);
    printf("\n");


    renormalize(A, B);

    printmat_complex(B, 3, 3);


}
