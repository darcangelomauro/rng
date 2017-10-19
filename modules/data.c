#define DATA_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "matop.h"
#include "global.h"


// PRINT SIMULATION DATA

void print_data(FILE* fdata)
{
    // print data for program
    fprintf(fdata, "%d %d %d %.15lf %.15lf %d %d %d\n", dim, nH, nL, SCALE, G, Ntherm, Nsw, GAP);

    // print data for human
    fprintf(fdata, "dim: %d\n", dim);
    fprintf(fdata, "nH: %d\n", nH);
    fprintf(fdata, "nL: %d\n", nL);
    fprintf(fdata, "SCALE: %lf\n", SCALE);
    fprintf(fdata, "G: %lf\n", G);
    fprintf(fdata, "Ntherm: %d\n", Ntherm);
    fprintf(fdata, "Nsw: %d\n", Nsw);
    fprintf(fdata, "GAP: %d\n", GAP);
}

void print_time(FILE* fdata, char* s)
{
    // format time
    time_t timer;
    char buffer[26];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    fprintf(fdata, "%s %s\n", s, buffer);
}


void print_thermalization(FILE* fobsS)
{
    // print action
    fprintf(fobsS, "%.15lf\n", S);
}

void print_thermalization_plus(FILE* fobsStr)
{
    fprintf(fobsStr, "%.15lf ", S);
    for(int i=0; i<nH; i++)
        fprintf(fobsStr, "%.15lf %.15lf ", tr[i], tr2[i]);
    fprintf(fobsStr, "\n");

}

void print_simulation(FILE* fobsS, FILE* fobsHL)
{
    // print action
    fprintf(fobsS, "%.15lf ", S);
    for(int i=0; i<nH; i++)
        fprintf(fobsS, "%.15lf %.15lf ", tr[i], tr2[i]);
    fprintf(fobsS, "\n");

    // print H matrices
    for(int l=0; l<nH; l++)
    {
        for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim; j++)
            {
                gsl_complex z = gsl_matrix_complex_get(H[l], i, j);
                fprintf(fobsHL, "%.15lf ", GSL_REAL(z)); 
                fprintf(fobsHL, "%.15lf ", GSL_IMAG(z));
            }
        }
        fprintf(fobsHL, "\n");
    }

    // print L matrices
    for(int l=0; l<nL; l++)
    {
        for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim; j++)
            {
                gsl_complex z = gsl_matrix_complex_get(L[l], i, j);
                fprintf(fobsHL, "%.15lf ", GSL_REAL(z)); 
                fprintf(fobsHL, "%.15lf ", GSL_IMAG(z));
            }
        }
        fprintf(fobsHL, "\n");
    }
}
