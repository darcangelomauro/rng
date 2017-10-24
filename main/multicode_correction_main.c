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

    for(int i=0; i<argc-1; i++)
    {
        // avoid overwriting previous analysis
        char* data = alloc_coded_filename("data", argv[i+1]);
        FILE* fdata = fopen(data, "r");
        char* simHL = alloc_coded_filename("simHL", argv[i+1]);
        FILE* fsimHL = fopen(simHL, "r");
        char* simS = alloc_coded_filename("simS", argv[i+1]);
        FILE* fsimS = fopen(simS, "r");
        char* simS_corr = alloc_coded_filename("simS2", argv[i+1]);
        FILE* fsimS_corr = fopen(simS_corr, "w");
        if(fdata == NULL || fsimHL == NULL || fsimS == NULL || fsimS_corr == NULL)
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


        int nmeas = Nsw_ / GAP_;
        for(int j=0; j<nmeas; j++)
        {
            double S = 0.;
            int r = fscanf(fsimS, "%lf", &S);
            for(int k=0; k<nH_; k++)
            {
                double temp;
                r += fscanf(fsimS, "%lf", &temp);
                r += fscanf(fsimS, "%lf", &temp);
            }

            if(r != 1+2*nH_)
            {
                printf("Error: failed to read action\n");
                exit(EXIT_FAILURE);
            }

            fprintf(fsimS_corr, "%.15lf ", S);
            
            
            
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
                fprintf(fsimS_corr, "%.15lf %.15lf ", trace_herm(H_[k]), trace_herm(H2_[k]));

            fprintf(fsimS_corr, "\n");
            

            
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

        
        // freeing memory
        free(data);
        fclose(fdata);
        free(simHL);
        fclose(fsimHL);
        free(simS);
        fclose(fsimS);
        free(simS_corr);
        fclose(fsimS_corr);
    }

}
