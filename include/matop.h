#ifndef MATOP_H
#define MATOP_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>

/* macro for disk memory mapping */
#define Mapmalloc(number, type, filename, fd)   \
    load_mmap((filename), &(fd), (number)*sizeof(type), 1)
#define Mapload(number, type, filename, fd)     \
    load_mmap((filename), &(fd), (number)*sizeof(type), 0)
#define Mapfree(number, type, fd, pointer)      \
    release_mmap((pointer), (number)*sizeof(type), (fd))


extern void printmat_complex(gsl_matrix_complex *m, int n, int k);
extern void fprintmat_complex(gsl_matrix_complex *m, int n, int k, FILE* f);
extern void dispherm(gsl_matrix_complex* m, int n);
extern void diag(gsl_matrix_complex* m, gsl_vector* eval, int n);
extern gsl_complex trace(gsl_matrix_complex* m);
extern gsl_complex trace_slow(gsl_matrix_complex* m, int n);
extern double trace_herm(gsl_matrix_complex* m);
extern void traceless(gsl_matrix_complex* m, int n);
extern void tensor(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c);
extern void tensor_refined(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c, gsl_complex v, int modeA, int modeB, int add);
extern void tensor_bruteforce(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c, int m, int n, int x, int y);
extern void tensor_sq(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c, int n);
extern void commutator(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void anticommutator(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void commutator_bruteforce(gsl_matrix_complex* m, gsl_matrix_complex* c, int n);
extern void anticommutator_bruteforce(gsl_matrix_complex* m, gsl_matrix_complex* c, int n);
extern void first_term(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void second_term(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void first_term_add(gsl_matrix_complex* m, gsl_matrix_complex* c, int mode);
extern void second_term_add(gsl_matrix_complex* m, gsl_matrix_complex* c, int mode);
extern void make_hermitian(gsl_matrix_complex* a);
#endif
