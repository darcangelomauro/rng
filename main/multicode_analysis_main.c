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

    FILE* fres = fopen("out.txt", "w");
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
        char* simS = alloc_coded_filename("simS", argv[i+1]);
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


        double* S;
        int nmeas = Nsw_ / GAP_;
        S = calloc(nmeas, sizeof(double));

        for(int j=0; j<nmeas; j++)
        {
            int r;
            r = fscanf(fsimS, "%lf", &S[j]);
            if(r==0)
            {
                printf("Error: not able to read action\n");
                exit(EXIT_FAILURE);
            }

            int r1 = 0;
            for(int k=0; k<2.*nH_; k++)
            {
                double temp;
                r1 += fscanf(fsimS, "%lf", &temp);
            }
            if(r1 != 2*nH_)
            {
                printf("Error: not enough data in %s\n", data);
                exit(EXIT_FAILURE);
            }
        }

        double mean = gsl_stats_mean(S, 1, nmeas);
        double var = gsl_stats_variance(S, 1, nmeas);

        fprintf(fres, "%lf %lf %lf\n", G_, mean, sqrt(var/(double)nmeas));
        free(S);
        free(data);
        fclose(fdata);
        free(simS);
        fclose(fsimS);
    }

    fclose(fres);

}
