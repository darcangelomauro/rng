#ifndef GLOBAL_H
#define GLOBAL_H

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>

// typical values:
//
// (1,0) GEOMETRY
// dim=5    SCALE=0.075
// dim=15   SCALE=0.013
//
// (2,0) GEOMETRY
// dim=5    SCALE=0.035
//
// (0,2) GEOMETRY
// dim=5     SCALE=0.035
// dim=10    SCALE=0.009
// dim=15    SCALE=0.0077
// dim=30    SCALE=0.0015
//
// (1,3) GEOMETRY
// dim=5, g=-3.65       SCALE=0.045 (acc. rate approx 0.25)
// dim=7, g=-3.70       SCALE=0.022 (acc. rate approx 0.25)


#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

// GLOBAL VARIABLES ********************************************************

#define GEOM04

EXTERN int dim;
#ifndef GAP
    #define GAP (100)
#endif
EXTERN int nH;
EXTERN int nL;
#ifndef nHL
    #define nHL (nH+nL) 
#endif
EXTERN double SCALE;        // metropolis scale factor
EXTERN double G;            // coupling constants
EXTERN double S;            // action
EXTERN int Nsw;             // number of simulation sweeps
EXTERN int Ntherm;          // number of thermalization sweeps

// H and L matrices
EXTERN gsl_matrix_complex** H;
EXTERN double* tr;     //this is for the traces of the H matrices
EXTERN double* tr2;    //this is fot the traces squared of the H matrices
EXTERN gsl_matrix_complex** L;          

// displacement matrix for metropolis move
// (just an auxiliary matrix but it is needed
EXTERN gsl_matrix_complex* M;

// dirac matrix for consistency check
EXTERN int dimG;
EXTERN int dimD;
EXTERN gsl_matrix_complex* DIRAC;
EXTERN gsl_matrix_complex** gammaH;
EXTERN gsl_matrix_complex** gammaL;          

#endif

