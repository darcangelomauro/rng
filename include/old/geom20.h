#ifndef GEOM20_H
#define GEOM20_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>

extern void geom_check20();
extern void init_gamma20();
extern double actionD2t20();
extern double actionD4D2t20();
extern double actionD4D2t20_realitycheck();
#endif
