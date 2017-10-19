#define UPDATE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_double.h>
#include "update.h"
#include "matop.h"
#include "fileop.h"
#include "math.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include "global.h"
#include "macros.h"
#include "data.h"
#include "statistics.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#define NAME(vec, idx, alt, check)  (idx == check ? (alt) : (vec[idx]))


// generates nxn hermitian matrix H (mode 0) or traceless hermitian matrix L (mode 1)
void generate_HL(gsl_matrix_complex* m, int mode, int n, gsl_rng* r)
{
    // generate random hermitian matrix
    for(int i=0; i<n; i++)
    {

        // off diagonal part
        for(int j=0; j<i; j++)
        {
            double x, y;

            // generate x and y uniformly between -1 and 1
            x = -1 + 2*gsl_rng_uniform(r);
            y = -1 + 2*gsl_rng_uniform(r);

            gsl_matrix_complex_set(m, i, j, gsl_complex_rect(x,y));
            gsl_matrix_complex_set(m, j, i, gsl_complex_rect(x,-y));
        }


        // diagonal part
        double z;
        z = -1 + 2*gsl_rng_uniform(r);
        gsl_matrix_complex_set(m, i, i, gsl_complex_rect(z,0.));

    }



    // turn traceless if mode 1
    if(mode)
    {
        //antihermitian we don't need, the i can be included at any time
        //gsl_matrix_complex_scale(m, gsl_complex_rect(0.,1.));

        //render traceless
        traceless(m, n);
    }


    
    // apply scale
    gsl_matrix_complex_scale(m, gsl_complex_rect(SCALE,0.));
    

    return;
}

void overwrite_init_file()
{
    FILE* finit = fopen("init.txt", "w");
    fprintf(finit, "%d\n", dim);
    fprintf(finit, "%d\n", nH);
    fprintf(finit, "%d\n", nL);
    fprintf(finit, "%lf\n", SCALE);
    fprintf(finit, "%lf\n", G);
    fprintf(finit, "%d\n", Ntherm);
    fprintf(finit, "%d\n", Nsw);
    fclose(finit);
}


void init_data()
{
    // read input data
    FILE* finit = fopen("init.txt", "r");
    if(finit == NULL)
    {
        printf("Error: input.txt not found\n");
        exit(EXIT_FAILURE);
    }

    int narg;
    // initialize matrix dimension
    narg = fscanf(finit, "%d", &dim);
    // initialize #H matrices
    narg += fscanf(finit, "%d", &nH);
    // initialize #L matrices
    narg += fscanf(finit, "%d", &nL);
    // initialize SCALE
    narg += fscanf(finit, "%lf", &SCALE);
    // initialize coupling constant
    narg += fscanf(finit, "%lf", &G);
    // initialize number of thermalization sweeps
    narg += fscanf(finit, "%d", &Ntherm);
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw);

    if(narg < 7)
    {
        printf("Error: not enough data in init.txt\n");
        exit(EXIT_FAILURE);
    }

    // initialize dimD
    if(nHL == 1)
        dimG = 1;
    else if(nHL<4)
        dimG = 2;
    else
        dimG = 4;
    dimD = dimG*dim*dim;
}

void init_data_analysis()
{
    // read input data
    FILE* finit = fopen("init_analysis.txt", "r");
    if(finit == NULL)
    {
        printf("Error: input_analysis.txt not found\n");
        exit(EXIT_FAILURE);
    }

    int narg;
    // initialize matrix dimension
    narg = fscanf(finit, "%d", &dim);
    // initialize #H matrices
    narg += fscanf(finit, "%d", &nH);
    // initialize #L matrices
    narg += fscanf(finit, "%d", &nL);
    // initialize SCALE
    narg += fscanf(finit, "%lf", &SCALE);
    // initialize coupling constant
    narg += fscanf(finit, "%lf", &G);
    // initialize number of thermalization sweeps
    narg += fscanf(finit, "%d", &Ntherm);
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw);

    if(narg < 7)
    {
        printf("Error: not enough data in init_analysis.txt\n");
        exit(EXIT_FAILURE);
    }

    // initialize dimD
    if(nHL == 1)
        dimG = 1;
    else if(nHL<5)
        dimG = 2;
    else
        dimG = 4;
    dimD = dimG*dim*dim;
}

//initializes H matrices to unity and L to zero
void init_minimal(void init_gamma())
{
    //initialize H matrices to unity
    H = malloc(nH*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        H[i] = gsl_matrix_complex_alloc(dim, dim);

    //initialize L matrices to zero
    L = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nL; i++)
        L[i] = gsl_matrix_complex_alloc(dim, dim);
    
    // allocate displacement matrix
    M = gsl_matrix_complex_alloc(dim, dim);
    
    // initialize gamma matrices and dirac
    gammaH = malloc(nH*sizeof(gsl_matrix_complex*));
    gammaL = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        gammaH[i] = gsl_matrix_complex_calloc(dimG, dimG);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
}
//initializes H matrices to unity and L to zero
void init_cold(void Sfunc(double*, int*), void init_gamma())
{
    //initialize H matrices to unity
    H = malloc(nH*sizeof(gsl_matrix_complex*));
    tr = malloc(nH*sizeof(double));
    tr2 = malloc(nH*sizeof(double));
    for(int i=0; i<nH; i++)
    {
        H[i] = gsl_matrix_complex_calloc(dim, dim);
        gsl_matrix_complex_set_identity(H[i]);
        tr[i] = dim;
        tr2[i] = dim;
    }

    //initialize L matrices to zero
    L = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nL; i++)
        L[i] = gsl_matrix_complex_calloc(dim, dim);
    

    // allocate displacement matrix
    M = gsl_matrix_complex_calloc(dim, dim);
    

    // initialize action
    double var[6];
    int control[6] = {1,0,0,0,0,0};
    Sfunc(var, control);
    S = var[0];


    // initialize gamma matrices and dirac
    gammaH = malloc(nH*sizeof(gsl_matrix_complex*));
    gammaL = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        gammaH[i] = gsl_matrix_complex_calloc(dimG, dimG);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
}

//initializes H matrices to unity and L to custom value
void init_custom(void Sfunc(double*, int*), void init_gamma(), char* filename)
{
    FILE* finitHL = fopen(filename, "r");
    //initialize H matrices
    H = malloc(nH*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
    {
        H[i] = gsl_matrix_complex_calloc(dim, dim);
        for(int j=0; j<dim*dim; j++)
        {
            double re, im;
            int read=0;
            read += fscanf(finitHL, "%lf", &re);
            read += fscanf(finitHL, "%lf", &im);
            if(read != 2)
            {
                printf("Error: custom initHL file not valid\n");
                exit(EXIT_FAILURE);
            }
            H[i] -> data[2*j] = re;
            H[i] -> data[2*j+1] = im;
        }
    }

    //initialize L matrices
    L = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nL; i++)
    {
        L[i] = gsl_matrix_complex_calloc(dim, dim);
        for(int j=0; j<dim*dim; j++)
        {
            double re, im;
            int read=0;
            read += fscanf(finitHL, "%lf", &re);
            read += fscanf(finitHL, "%lf", &im);
            if(read != 2)
            {
                printf("Error: custom initHL file not valid\n");
                exit(EXIT_FAILURE);
            }
            L[i] -> data[2*j] = re;
            L[i] -> data[2*j+1] = im;
        }
    }

    fclose(finitHL);
    

    // allocate displacement matrix
    M = gsl_matrix_complex_calloc(dim, dim);
    

    // initialize action and traces
    double var[6];
    int control[6] = {1,0,1,1,1,1};
    Sfunc(var, control);

    S = var[0];
    for(int i=0; i<nH; i++)
    {
        tr[i] = var[nH+2];
        tr2[i] = var[nH+3];
    }


    // initialize gamma matrices and dirac
    gammaH = malloc(nH*sizeof(gsl_matrix_complex*));
    gammaL = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        gammaH[i] = gsl_matrix_complex_calloc(dimG, dimG);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
}

//initializes H and L with random matrices
void init_hot(void Sfunc(double*, int*), void init_gamma(), gsl_rng* r)
{
    //initialize H matrices randomly
    H = malloc(nH*sizeof(gsl_matrix_complex*));
    tr = malloc(nH*sizeof(double));
    tr2 = malloc(nH*sizeof(double));
    for(int i=0; i<nH; i++)
    {
        H[i] = gsl_matrix_complex_calloc(dim, dim);
        generate_HL(H[i], 0, dim, r);
    }
    
    //initialize L matrices randomly
    L = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nL; i++)
    {
        L[i] = gsl_matrix_complex_calloc(dim, dim);
        generate_HL(L[i], 1, dim, r);
    }


    // allocate displacement matrix
    M = gsl_matrix_complex_calloc(dim, dim);
    

    // initialize action and traces
    double var[6];
    int control[6] = {1,0,1,1,1,1};
    Sfunc(var, control);

    S = var[0];
    for(int i=0; i<nH; i++)
    {
        tr[i] = var[nH+2];
        tr2[i] = var[nH+3];
    }

    
    // initialize gamma matrices
    gammaH = malloc(nH*sizeof(gsl_matrix_complex*));
    gammaL = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        gammaH[i] = gsl_matrix_complex_calloc(dimG, dimG);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    // alloc dirac
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
}


void simulation_free()
{
    // free H matrices
    for(int i=0; i<nH; i++)
        gsl_matrix_complex_free(H[i]);
    free(H);
    free(tr);
    free(tr2);

    // free L matrices
    for(int i=0; i<nL; i++)
        gsl_matrix_complex_free(L[i]);
    free(L);

    // free displacement matrix
    gsl_matrix_complex_free(M);
    
    // free gamma matrices and dirac
    for(int i=0; i<nH; i++)
        gsl_matrix_complex_free(gammaH[i]);
    free(gammaH);
    for(int i=0; i<nL; i++)
        gsl_matrix_complex_free(gammaL[i]);
    free(gammaL);
    gsl_matrix_complex_free(DIRAC);
}



// returns 1 if ith sweep corresponds to
// a measurement (step is the distance between
// measurements)
int measurement(int i, int step)
{
    return !(i%step);
}


int move(void Sfunc(double*, int*), int mode, gsl_rng* r)
{
    // check is what will be returned
    int check = 0;

    // MONTECARLO MOVE PROPOSAL
    double S1[14];
    int control[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0};
    if(mode) control[1] = 1;

    
    // decide what matrix gets updated
    int uM = (int)(nHL*gsl_rng_uniform(r));
    while(uM == nHL)
        uM = (int)(nHL*gsl_rng_uniform(r));

    
    // buffer is needed to temporarily update the H or L matrices
    // but it never allocates new memory, so no need to free at the end
    gsl_matrix_complex* buffer;
    

    if(uM < nH)
    {
        // generate displacement dH
        generate_HL(M, 0, dim, r);    //M is now the proposed displacement from H[uM]
        control[uM+2] = 1;
        control[uM+3] = 1;

        // compute displaced matrix
        gsl_matrix_complex_add(M, H[uM]);             //now M = H[uM] + displacement
        buffer = H[uM];
        H[uM] = M;
        Sfunc(S1, control);
        H[uM] = buffer;
    }

    else if(uM<nHL)
    {
        // generate displacement dL
        generate_HL(M, 1, dim, r);    //M is now the proposed displacement from L[uM-nH]

        // compute displaced matrix
        gsl_matrix_complex_add(M, L[uM-nH]);                  //now L[uM-nH] = L[uM-nH] + displacement
        buffer = L[uM-nH];
        L[uM-nH] = M;
        Sfunc(S1, control);
        L[uM-nH] = buffer;
    }

    else
    {
        printf("error: index out of bounds\n");
        exit(EXIT_FAILURE);
    }


    // now we evaluate the move
    if(S1[0] < S)
    {
        // update H or L
        if(uM < nH)
        {
            gsl_matrix_complex_memcpy(H[uM], M);
            tr[uM] = S1[uM+2];
            tr2[uM] = S1[uM+3];
        }
        else if(uM<nHL)
            gsl_matrix_complex_memcpy(L[uM-nH], M);
        else
        {
            printf("error: index out of bounds\n");
            exit(EXIT_FAILURE);
        }

        // update action
        S = S1[0];

        // move accepted
        check = 1;
    }
    else
    {
        double e = exp(S-S1[0]);
        double p = gsl_rng_uniform(r);

        if(e>p)
        {
            // update H or L
            if(uM < nH)
            {
                gsl_matrix_complex_memcpy(H[uM], M);
                tr[uM] = S1[uM+2];
                tr2[uM] = S1[uM+3];
            }
            else if(uM<nHL)
                gsl_matrix_complex_memcpy(L[uM-nH], M);
            else
            {   
                printf("error: index out of bounds\n");
                exit(EXIT_FAILURE);
            }   

            // update action
            S = S1[0];

            // move accepted
            check = 1;
        }
    }

    return check;

}

// a sweep is a collection of nHL proposed moves
// it returns the acceptance rate of that sweep
double sweep(void Sfunc(double*, int*), int mode, gsl_rng* r)
{
    double sum = 0.;
    for(int i=0; i<nHL; i++)
        sum += move(Sfunc, mode, r);


    return sum/(double)nHL;
}


// a routine to automatically set the SCALE factor to give acceptance rate
// between minTarget and maxTarget (very rudimental and ugly, seems to work reasonably well)
void SCALE_autotune(double minTarget, double maxTarget, void Sfunc(double*, int*), gsl_rng* r)
{
    int n1 = 50;
    int n2 = 50;
    int m = 50;
    double limit = 1e-5;
    double ar = 0.;

    // initialize ar
    for(int i=0; i<n1; i++)
        ar += sweep(Sfunc, 0, r);
    ar /= (double)n1;

    int count = 0;

    while((ar < minTarget || ar > maxTarget) && count < 2 )
    {
        count++;

        // tune delta
        double temp = SCALE;
        int order;
        for(order=0; order<1000; order++)
        {
            if((int)temp != 0)
                break;
            else
                temp *= 10.;
        }
        double delta = pow(10., -order);

        // tune SCALE
        while(delta > limit)
        {
            for(int j=0; j<m; j++)
            {
                if(ar > maxTarget)
                    SCALE += delta;
                else if(ar < minTarget)
                {
                    if((SCALE-delta) > 0)
                        SCALE -= delta;
                    else
                        break;
                }
                ar = 0.;
                for(int k=0; k<n2; k++)
                    ar += sweep(Sfunc, 0, r);
                ar /= (double)n2;
            }
            delta /= 10.;
        }
    }
}


// complete simulation routine
// first, it prints on file the simulation data
// the thermalization part outputs two files "thermalizationX.txt" (where X is 1 or 2) with the action value
// the simulation part outputs a file "simulation.txt" with the action and the H and L matrices
// returns acceptance rate
double simulation(void Sfunc(double*, int*), int mode, void init_gamma(), gsl_rng* r)
{
    init_data();
    
    GEOM_CHECK();
    
    init_cold(Sfunc, init_gamma);
    SCALE_autotune(0.21, 0.3, Sfunc, r);
    printf("Auto tuned SCALE: %lf\n", SCALE);

    // generate unique filename
    char* code = generate_code(5, r);
    printf("%s\n", code);
    char* name_data = alloc_coded_filename("data", code);
    char* name_therm1 = alloc_coded_filename("therm1", code);
    char* name_therm2 = alloc_coded_filename("therm2", code);
    char* name_simS = alloc_coded_filename("simS", code);
    char* name_simHL = alloc_coded_filename("simHL", code);

    FILE* fdata = fopen(name_data, "w");
    FILE* ftherm1 = fopen(name_therm1, "w");
    FILE* ftherm2 = fopen(name_therm2, "w");
    FILE* fsimS = fopen(name_simS, "w");
    FILE* fsimHL = fopen(name_simHL, "w");

    free(code);
    free(name_data);
    free(name_therm1);
    free(name_therm2);
    free(name_simS);
    free(name_simHL);

    // print simulation data
    print_data(fdata);

    // thermalization1
    init_cold(Sfunc, init_gamma);
    print_time(fdata, "start therm1:");
    for(int i=0; i<Ntherm; i++)
    {
        sweep(Sfunc, mode, r);
        print_thermalization_plus(ftherm1);
    }
    print_time(fdata, "end therm1:");
    fclose(ftherm1);
    simulation_free();
    
    // thermalization2
    init_hot(Sfunc, init_gamma, r);
    print_time(fdata, "start therm2:");
    for(int i=0; i<Ntherm; i++)
    {
        sweep(Sfunc, mode, r);
        print_thermalization_plus(ftherm2);
    }
    print_time(fdata, "end therm2:");
    fclose(ftherm2);

    
    // simulation (starts from second thermalization)
    double ar = 0.;
    print_time(fdata, "start simulation:");
    for(int i=0; i<Nsw; i++)
    {
        ar += sweep(Sfunc, mode, r);
        if(measurement(i, GAP))
            print_simulation(fsimS, fsimHL);
    }
    print_time(fdata, "end simulation:");
    fprintf(fdata, "acceptance rate: %lf", ar/(double)Nsw);

    fclose(fsimS);
    fclose(fsimHL);
    fclose(fdata);
    simulation_free();

    return ar/(double)Nsw;
}

void analysis(char* code, int* control, void init_gamma())
{
    /*
    char* data;
    char* simM;
    int n = strlen(code);
    data = malloc((n+10)*sizeof(char));
    simM = malloc((n+10)*sizeof(char));
    alloc_input_filename(data, simM, code);
    */
    char* data = alloc_coded_filename("data", code);
    char* simM = alloc_coded_filename("simM", code);

    FILE* fdata = fopen(data, "r");
    if(fdata == NULL)
    {
        printf("Error: %s not found\n", data);
        exit(EXIT_FAILURE);
    }
    FILE* fsimM = fopen(simM, "r");
    if(fsimM == NULL)
    {
        printf("Error: %s not found\n", simM);
        exit(EXIT_FAILURE);
    }
    free(data);
    free(simM);
    FILE* finit = fopen("init_analysis.txt", "w");

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

    fclose(fdata);

    fprintf(finit, "%d\n", dim_);
    fprintf(finit, "%d\n", nH_);
    fprintf(finit, "%d\n", nL_);
    fprintf(finit, "%.15lf\n", SCALE_);
    fprintf(finit, "%.15lf\n", G_);
    fprintf(finit, "%d\n", Ntherm_);
    fprintf(finit, "%d\n", Nsw_);
    fclose(finit);
    printf("init_analysis.txt written\n");

    init_data_analysis();
    init_minimal(init_gamma);
    
    GEOM_CHECK();
    
    int N = Nsw/GAP_;
    /*
    char* evlM;
    char* evlD;
    char* asys;
    evlM = malloc((n+10)*sizeof(char));
    evlD = malloc((n+10)*sizeof(char));
    asys = malloc((n+10)*sizeof(char));
    alloc_input_filename2(evlM, evlD, asys, code);
    */
    char* evalM = alloc_coded_filename("evalM", code);
    char* evalD = alloc_coded_filename("evalD", code);
    char* asys = alloc_coded_filename("asys_time", code);
    FILE* fevalM = fopen(evalM, "w");
    FILE* fevalD = fopen(evalD, "w");
    FILE* fasys = fopen(asys, "w");
    free(evalM);
    free(evalD);
    free(asys);
    print_time(fasys, "start analisys:");
    for(int k=0; k<N; k++)
    {
        if(!(k%(N/10)))
            printf("processing sweep #: %d\n", k);

        // read H matrices
        for(int l=0; l<nH; l++)
        {
            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    double re, im;
                    int nmatel = 0;
                    nmatel += fscanf(fsimM, "%lf", &re);
                    nmatel += fscanf(fsimM, "%lf", &im);
                    if(nmatel != 2)
                    {
                        printf("Error: not enough matrices\n");
                        exit(EXIT_FAILURE);
                    }
                    gsl_matrix_complex_set(H[l], i, j, gsl_complex_rect(re, im));
                }
            }
        }
        // read L matrices
        for(int l=0; l<nL; l++)
        {
            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    double re, im;
                    int nmatel = 0;
                    nmatel += fscanf(fsimM, "%lf", &re);
                    nmatel += fscanf(fsimM, "%lf", &im);
                    if(nmatel != 2)
                    {
                        printf("Error: not enough matrices\n");
                        printf("k = %d\n", k);
                        exit(EXIT_FAILURE);
                    }
                    
                    gsl_matrix_complex_set(L[l], i, j, gsl_complex_rect(re, im));
                }
            }
        }


        // now the H and L matrices are set, and we can compute observables

        if(control[0])
        {
            gsl_vector* vevalH = gsl_vector_alloc(dim);
            for(int i=0; i<nH; i++)
            {
                diag(H[i], vevalH, dim);
                for(int j=0; j<dim; j++)
                    fprintf(fevalM, "%.15lf ", gsl_vector_get(vevalH, j));
                fprintf(fevalM, "\n");
            }
            gsl_vector_free(vevalH);
        }
        
        if(control[1])
        {
            gsl_vector* vevalL = gsl_vector_alloc(dim);
            for(int i=0; i<nL; i++)
            {
                diag(L[i], vevalL, dim);
                for(int j=0; j<dim; j++)
                    fprintf(fevalM, "%.15lf ", gsl_vector_get(vevalL, j));
                fprintf(fevalM, "\n");
            }
            gsl_vector_free(vevalL);
        }
        
        if(control[2])
        {
            gsl_vector* vevalD = gsl_vector_alloc(dimD);
            build_dirac();
            diag(DIRAC, vevalD, dimD);
            for(int i=0; i<dimD; i++)
                fprintf(fevalD, "%.15lf ", gsl_vector_get(vevalD, i));
            fprintf(fevalD, "\n");

            gsl_vector_free(vevalD);
        }
    }
    print_time(fasys, "end analisys:");
    fclose(fsimM);
    fclose(fevalM);
    fclose(fevalD);
    fclose(fasys);
    simulation_free();
}

void hermitization()
{
    for(int i=0; i<nH; i++)
        make_hermitian(H[i]);
    for(int i=0; i<nL; i++)
        make_hermitian(L[i]);
}

