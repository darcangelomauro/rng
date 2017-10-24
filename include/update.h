#ifndef UPDATE_H
#define UPDATE_H

#include "global.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_rng.h>


extern void generate_HL(gsl_matrix_complex* m, int mode, int n, gsl_rng* r);
extern void overwrite_init_file();
extern void init_data();
extern void init_data_analysis();
extern void init_minimal(void init_gamma());
extern void init_cold(void Sfunc(double*, int*), void init_gamma());
extern void init_custom(void Sfunc(double*, int*), void init_gamma(), char* filename);
extern void init_hot(void Sfunc(double*, int*), void init_gamma(), gsl_rng* r);
extern void simulation_free();
extern int measurement(int i, int step);
extern int move(void Sfunc(double*, int*), int mode, gsl_rng* r);
extern double sweep(void Sfunc(double*, int*), int mode, gsl_rng* r);
extern void SCALE_autotune(double minTarget, double maxTarget, void Sfunc(double*, int*), gsl_rng* r);
extern char* simulation(void Sfunc(double*, int*), int mode, void init_gamma(), gsl_rng* r);
extern void analysis(char* code, int* control, void init_gamma());
extern void hermitization();
#endif
