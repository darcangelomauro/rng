#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include <gsl/gsl_histogram.h>
#include "update.h"
#include "global.h"
#include "macros.h"

#define min -1.5
#define max 1.5
#define nbin 101
#define N 10000
#define dimbin ((max-min)/(1.*nbin))


int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("gimme the code\n");
        exit(EXIT_FAILURE);
    }
    
    char* evalM = alloc_coded_filename("evalM", argv[1]);
    char* evalD = alloc_coded_filename("evalD", argv[1]);
    FILE* fevalM = fopen(evalM, "r");
    FILE* fevalD = fopen(evalD, "r");
    free(evalM);
    free(evalD);

    gsl_histogram* histoM = gsl_histogram_alloc(nbin);
    gsl_histogram_set_ranges_uniform(histoM, min, max);
    gsl_histogram* histoD = gsl_histogram_alloc(nbin);
    gsl_histogram_set_ranges_uniform(histoD, min, max);

    //while(0)
    for(int i=0; i<N; i++)
    {
        double val;
        int read = fscanf(fevalM, "%lf", &val);
        if(read)
            gsl_histogram_increment(histoM, val);
        else
            break;
    }
    fclose(fevalM);

    while(0)
    //for(int i=0; i<N; i++)
    {
        double val;
        int read = fscanf(fevalD, "%lf", &val);
        if(read)
            gsl_histogram_increment(histoD, val);
    }
    fclose(fevalD);
    
    gsl_histogram_scale(histoM, 1./(N*dimbin));
    gsl_histogram_scale(histoD, 1./(N*dimbin));


    char* chistoM = alloc_coded_filename("histoM", argv[1]);
    char* chistoD = alloc_coded_filename("histoD", argv[1]);
    FILE* foutM = fopen(chistoM, "w");
    FILE* foutD = fopen(chistoD, "w");
    free(chistoM);
    free(chistoD);
    gsl_histogram_fprintf(foutM, histoM, "%.15lf", "%.15lf");
    gsl_histogram_fprintf(foutD, histoD, "%.15lf", "%.15lf");
    fclose(foutM);
    fclose(foutD);
    gsl_histogram_free(histoM);
    gsl_histogram_free(histoD);
}
