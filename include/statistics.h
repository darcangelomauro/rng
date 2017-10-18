#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdlib.h>
#include <stdio.h>

//extern void alloc_input_filename(char* data, char* simM, char* code);
//extern void alloc_input_filename2(char* evlM, char* evlD, char* asys, char* code);
extern void read_input(FILE* fsimM, int dim_s, int nH_s, int nL_s, int n, FILE* fout);
extern double* binned_vector(double* v, int n, int dimbin);
extern int binned_size(int n, int dimbin);
#endif
