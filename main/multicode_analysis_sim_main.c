#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <fileop.h>
#include <matop.h>
#include <global.h>

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        printf("gimme some code\n");
        exit(EXIT_FAILURE);
    }

    FILE* fres = fopen("out_sim.txt", "w");
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
        char* simHL = alloc_coded_filename("simHL", argv[i+1]);
        FILE* fsimHL = fopen(simHL, "r");
        if(fdata == NULL || fsimHL == NULL)
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
            gsl_matrix_complex** H_ = malloc(nH_*sizeof(gsl_matrix_complex*));
            gsl_matrix_complex** H2_ = malloc(nH_*sizeof(gsl_matrix_complex*));
            gsl_matrix_complex** L_ = malloc(nL_*sizeof(gsl_matrix_complex*));
            for(int k=0; k<nH_; k++)
                H_[k] = gsl_matrix_complex_calloc(dim_, dim_);
            for(int k=0; k<nL_; k++)
                L_[k] = gsl_matrix_complex_calloc(dim_, dim_);
            
            
            int r1 = 0;
            for(int k=0; k<nH_; k++)
            {
                for(int l=0; l<dim_; l++)
                {
                    for(int m=0; m<dim_; m++)
                    {
                        double re = 0.;
                        double im = 0.;
                        r1 += fscanf(fsimHL, "%lf", &re);
                        r1 += fscanf(fsimHL, "%lf", &im);
                        gsl_matrix_complex_set(H_[k], l, m, gsl_complex_rect(re, im));
                    }
                }
            }
            for(int k=0; k<nL_; k++)
            {
                for(int l=0; l<dim_; l++)
                {
                    for(int m=0; m<dim_; m++)
                    {
                        double re = 0.;
                        double im = 0.;
                        r1 += fscanf(fsimHL, "%lf", &re);
                        r1 += fscanf(fsimHL, "%lf", &im);
                        gsl_matrix_complex_set(L_[k], l, m, gsl_complex_rect(re, im));
                    }
                }
            }
            if(r1 != (nH_+nL_)*dim_*dim_*2)
            {
                printf("Error: not enough data in %s\n", simHL);
                exit(EXIT_FAILURE);
            }
            
            for(int k=0; k<nH_; k++)
            {
                H2_[k] = gsl_matrix_complex_calloc(dim_, dim_);
                gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, H_[k], H_[k], GSL_COMPLEX_ZERO, H2_[k]);
            }


            for(int k=0; k<nH_; k++)
            {
                tr[k][j] = trace_herm(H_[k]);
                tr2[k][j] = trace_herm(H2_[k]);
            }


            
            // freeing memory
            for(int k=0; k<nH_; k++)
                gsl_matrix_complex_free(H_[k]);
            for(int k=0; k<nH_; k++)
                gsl_matrix_complex_free(H2_[k]);
            for(int k=0; k<nL_; k++)
                gsl_matrix_complex_free(L_[k]);
            free(H_);
            free(H2_);
            free(L_);
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


        
        // freeing memory
        for(int j=0; j<nH_; j++)
        {
            free(tr[j]);
            free(tr2[j]);
        }
        free(tr);
        free(tr2);
        free(data);
        fclose(fdata);
        free(simHL);
        fclose(fsimHL);
    }

    fclose(fres);

}
