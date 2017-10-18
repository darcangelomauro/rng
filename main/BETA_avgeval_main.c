#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include <gsl/gsl_histogram.h>
#include "update.h"
#include "global.h"
#include "macros.h"

#define mindimbin 10
#define maxdimbin 1000
#define ntrials (maxdimbin-mindimbin)

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        printf("gimme the code and the D or M option\n");
        exit(EXIT_FAILURE);
    }
    if(argc < 3)
    {
        printf("gimme the D or M option\n");
        exit(EXIT_FAILURE);
    }

    int check;
    if(argv[2] == "D")
        check = 1;
    else if(argv[2] == "M")
        check = 0;
    else
    {
        printf("Option not valid. write D if the eigenvalues are Dirac or M if they are H or L\n");
        exit(EXIT_FAILURE);
    }
    
    char* eval;
    if(check)
        eval = alloc_coded_filename("evalD", argv[1]);
    else
        eval = alloc_coded_filename("evalM", argv[1]);

    char* data = alloc_coded_filename("data", argv[1]);
    FILE* feval = fopen(eval, "r");
    FILE* fdata = fopen(data, "r");
    free(eval);

    int narg=0;
    int dim_, nH_, nL_, Ntherm_, Nsw_, GAP_, dimG_, nHL_;
    double SCALE_, G_;
    // initialize matrix dimension
    narg += fscanf(fdata, "%d", &dim_);
    // initialize #H matrices
    narg += fscanf(fdata, "%d", &nH_);
    // initialize #L matrices
    narg += fscanf(fdata, "%d", &nL_);
    // initialize SCALE
    narg += fscanf(fdata, "%lf", &SCALE_);
    // initialize coupling constant
    narg += fscanf(fdata, "%lf", &G_);
    // initialize number of thermalization sweeps
    narg += fscanf(fdata, "%d", &Ntherm_);
    // initialize number of simulation sweeps
    narg += fscanf(fdata, "%d", &Nsw_);
    // initialize GAP
    narg += fscanf(fdata, "%d", &GAP_);

    if(narg < 8)
    {
        printf("Error: not enough data in %s\n", data);
        exit(EXIT_FAILURE);
    }

    fclose(fdata);
    free(data);

    nHL_ = nH_+nL_;
    if(nHL_ == 1)
        dimG_ = 1;
    else if(nHL_<5)
        dimG_ = 2;
    else
        dimG_ = 4;

    // number of eigenvalues
    int neval;
    if(check)
        neval = dim_*dim_*dimG_;
    else
        neval = dim_;

    // number of measurements
    int nmeas = Nsw_/GAP_;

    double** values = malloc(neval*sizeof(double*));
    for(int i=0; i<neval; i++)
        values[i] = malloc(nmeas*sizeof(double));

    // now valuesD and valuesM are a (neval)x(nmeas) matrix
    // that has to be initialized with the measurements
    for(int i=0; i<nmeas; i++)
    {
        for(int j=0; j<neval; j++)
        {
            int read;
            read = fscanf(feval, "%lf", &values[j][i]);
            if(read != 1)
            {
                printf("Error, failed to read evals\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    char* out = alloc_coded_filename("avgevalD", argv[1]);
    FILE* fout = fopen(out, "w");

    double** variance = malloc(ntrials*sizeof(double*));
    for(int i=0; i<ntrials; i++)
        variance[i] = malloc(neval*sizeof(double));
    for(int k=mindimbin; k<maxdimbin; k++)
    {
        for(int i=0; i<neval; i++)
        {
            // i is the eigenvalue index
            // k is the dimbin index
            double* binned = binned vector(values[i], nmeas, k); 










