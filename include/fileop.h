#ifndef FILEOP_H
#define FILEOP_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

extern void build_filename(char* name, char* directory, char* extension, double num);
extern char* generate_code(int n, gsl_rng* r);
extern char* alloc_coded_filename(char* suffix, char* code);
extern void avgfile(int row, int col, char *input, char *output);
extern void DfromH_t10(int row, int col, char *input, char *output);
extern void histo(char* inname, char* outname, int Nbin, int min, int max, int row, int col);
#endif
