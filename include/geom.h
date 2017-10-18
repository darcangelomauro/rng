#ifndef GEOM_H
#define GEOM_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>

#ifndef GEOM10_C
extern void geom_check10();
extern void init_gamma10();
extern void actionD2t10(double* vecS, int* control);
extern void actionD4D2t10(double* vecS, int* control);
#endif

#ifndef GEOM01_C
extern void geom_check01();
extern void init_gamma01();
extern void actionD2t01(double* vecS, int* control);
extern void actionD4D2t01(double* vecS, int* control);
#endif

#ifndef GEOM20_C
extern void geom_check20();
extern void init_gamma20();
extern void actionD2t20(double* vecS, int* control);
extern void actionD4D2t20(double* vecS, int* control);
#endif

#ifndef GEOM02_C
extern void geom_check02();
extern void init_gamma02();
extern void actionD2t02(double* vecS, int* control);
extern void actionD4D2t02(double* vecS, int* control);
#endif

#ifndef GEOM11_C
extern void geom_check11();
extern void init_gamma11();
extern void actionD2t11(double* vecS, int* control);
extern void actionD4D2t11(double* vecS, int* control);
#endif

#ifndef GEOM03_C
extern void geom_check03();
extern void init_gamma03();
extern void actionD2t03(double* vecS, int* control);
extern void actionD4D2t03(double* vecS, int* control);
#endif

#ifndef GEOM13_C
extern void geom_check13();
extern void init_gamma13();
extern void actionD2t13(double* vecS, int* control);
extern void actionD4D2t13(double* vecS, int* control);
#endif


#endif
