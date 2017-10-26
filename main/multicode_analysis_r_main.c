#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <fileop.h>
#include <global.h>

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        printf("gimme some code\n");
        exit(EXIT_FAILURE);
    }

    FILE* fres = fopen("out_r.txt", "w");
    if(fres == NULL)
    {
        printf("Error: unable to open output file\n");
        exit(EXIT_FAILURE);
    }

    for(int i=0; i<argc-1; i++)
    {
        // avoid overwriting previous analysis
        char* data = alloc_coded_filename("data", argv[i+1]);
        FILE* fdata = fopen(data, "r");
        char* simS = alloc_coded_filename("simS_r", argv[i+1]);
        FILE* fsimS = fopen(simS, "r");
        if(fdata == NULL || fsimS == NULL)
        {
            printf("Error: unable to read input file %s\n", argv[i+1]);
            exit(EXIT_FAILURE);
        }


        int narg=0;
        int dim_, nH_, nL_, Ntherm_, Nsw_, GAP_;
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
        dim_ /= 2;


        double** tr;
        double** tr2;
        int nmeas = Nsw_ / GAP_;
        tr = malloc(nH_*sizeof(double*));
        tr2 = malloc(nH_*sizeof(double*));
        for(int j=0; j<nH_; j++)
        {
            tr[j] = calloc(nmeas, sizeof(double));
            tr2[j] = calloc(nmeas, sizeof(double));
        }

        for(int j=0; j<nmeas; j++)
        {
            int r1 = 0;
            for(int k=0; k<nH_; k++)
            {
                r1 += fscanf(fsimS, "%lf", &tr[k][j]);
                r1 += fscanf(fsimS, "%lf", &tr2[k][j]);
            }
            if(r1 != 2*nH_)
            {
                printf("Error: not enough data in %s\n", simS);
                exit(EXIT_FAILURE);
            }
        }

        fprintf(fres, "%lf ", G_);
        for(int j=0; j<nH_; j++)
        {
            double meantr = gsl_stats_mean(tr[j], 1, nmeas);
            double vartr = gsl_stats_variance(tr[j], 1, nmeas);
            double meantr2 = gsl_stats_mean(tr2[j], 1, nmeas);
            double vartr2 = gsl_stats_variance(tr2[j], 1, nmeas);
            fprintf(fres, "%lf %lf %lf %lf ", meantr, sqrt(vartr/(double)nmeas), meantr2, sqrt(vartr2/(double)nmeas));
        }
        fprintf(fres, "\n");

        for(int j=0; j<nH_; j++)
        {
            free(tr[j]);
            free(tr2[j]);
        }
        free(tr);
        free(tr2);
        free(data);
        fclose(fdata);
        free(simS);
        fclose(fsimS);
    }

    fclose(fres);

}
