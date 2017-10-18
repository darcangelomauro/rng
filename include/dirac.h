#ifndef DIRAC_H
#define DIRAC_H

#include <gsl/gsl_complex.h>

extern void build_dirac();
//extern void build_dirac_fast();
extern gsl_complex actionD4D2_bruteforce();
extern gsl_complex actionD2_bruteforce();
#endif
