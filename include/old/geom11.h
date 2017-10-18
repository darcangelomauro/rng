#ifndef GEOM11_H
#define GEOM11_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>


extern void geom_check11();
extern void init_gamma11();
extern double actionD2t11();
extern double actionD4D2t11();
extern double actionD4D2t11_realitycheck();
#endif
