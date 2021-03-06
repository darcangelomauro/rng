#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include "fileop.h"
#include "statistics.h"
#include "global.h"

#define DIMBIN 1000

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        printf("gimme some code\n");
        exit(EXIT_FAILURE);
    }

    FILE* fres = fopen("out_pq20_dim10_r_4e7.txt", "w");
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
        //dim_ /= 2;


        double* obs;
        double* bin_obs;
        double** tr;
        double** bin_tr;
        double** tr2;
        double** bin_tr2;
        int nmeas = Nsw_ / GAP_;
        int dimbin = DIMBIN;
        int nbin = binned_size(nmeas, dimbin);
        obs = calloc(nmeas, sizeof(double));
        tr = malloc(nH_*sizeof(double*));
        tr2 = malloc(nH_*sizeof(double*));
        bin_tr = malloc(nH_*sizeof(double*));
        bin_tr2 = malloc(nH_*sizeof(double*));
        for(int j=0; j<nH_; j++)
        {
            tr[j] = calloc(nmeas, sizeof(double));
            tr2[j] = calloc(nmeas, sizeof(double));
        }

        for(int j=0; j<nmeas; j++)
        {
            int r1 = 0;
            double norm = 0;
            for(int k=0; k<nH_; k++)
            {
                r1 += fscanf(fsimS, "%lf", &tr[k][j]);
                r1 += fscanf(fsimS, "%lf", &tr2[k][j]);
                obs[j] += tr[k][j]*tr[k][j];
                norm += tr2[k][j];
            }
            obs[j] /= dim_*norm;
            if(r1 != 2*nH_)
            {
                printf("Error: not enough data in %s\n", simS);
                exit(EXIT_FAILURE);
            }
        }

        bin_obs = binned_vector(obs, nmeas, dimbin);
        for(int j=0; j<nH_; j++)
        {
            bin_tr[j] = binned_vector(tr[j], nmeas, dimbin);
            bin_tr2[j] = binned_vector(tr2[j], nmeas, dimbin);
        }
        
        double meanobs = gsl_stats_mean(bin_obs, 1, nbin);
        double varobs = gsl_stats_variance(bin_obs, 1, nbin);
        fprintf(fres, "%lf %.15lf %.15lf ", G_, meanobs, sqrt(varobs/(double)nbin));
        for(int j=0; j<nH_; j++)
        {
            double meantr = gsl_stats_mean(bin_tr[j], 1, nbin);
            double vartr = gsl_stats_variance(bin_tr[j], 1, nbin);
            double meantr2 = gsl_stats_mean(bin_tr2[j], 1, nbin);
            double vartr2 = gsl_stats_variance(bin_tr2[j], 1, nbin);
            fprintf(fres, "%.15lf %.15lf %.15lf %.15lf ", meantr, sqrt(vartr/(double)nbin), meantr2, sqrt(vartr2/(double)nbin));
        }
        fprintf(fres, "\n");

        free(obs);
        free(bin_obs);
        for(int j=0; j<nH_; j++)
        {
            free(tr[j]);
            free(tr2[j]);
            free(bin_tr[j]);
            free(bin_tr2[j]);
        }
        free(tr);
        free(tr2);
        free(bin_tr);
        free(bin_tr2);
        free(data);
        fclose(fdata);
        free(simS);
        fclose(fsimS);
    }

    fclose(fres);

}
