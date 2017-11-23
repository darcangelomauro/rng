#ifndef MACROS_H
#define MACROS_H

// this beautiful header makes life easier and writing shorter.
// just #define the geometry type in the main program, and
// all the relevant #include will be included.
// moreover, when calling of passing functions for computing the
// action, the gamma matrices or the dirac operator, you no longer need to
// speficy the type of geometry.
// eg1: #define GEOM20 in your main, and then instead of writing dirac20(), just write C_dirac
// eg1: #define GEOM20 in your main, and then instead of writing dirac20, just write P_dirac
// from the examples should be clear that when you call a
// function you should use the prefix C_, while when you
// pass the function as an argument you should use the
// prefix P_

#include "matop.h"
#include "fileop.h"
#include "update.h"
#include "dirac.h"
#include "data.h"
#include "geom.h"

#define P_dirac build_dirac
#define P_actionD2_b actionD2_bruteforce;
#define P_actionD4D2_b actionD4D2_bruteforce;
#define C_dirac build_dirac()
#define C_actionD2_b actionD2_bruteforce();
#define C_actionD4D2_b actionD4D2_bruteforce();
//#define P_dirac_b build_dirac_bruteforce
//#define C_dirac_b build_dirac_bruteforce()


#ifdef GEOM20
#define cliff_p 2
#define cliff_q 0

#define P_actionD2 actionD2t20
#define P_actionD4D2 actionD4D2t20
#define P_actionD4D2_r actionD4D2t20_debug
#define C_actionD2 actionD2t20()
#define C_actionD4D2 actionD4D2t20()
#define C_actionD4D2_r actionD4D2t20_debug()

#define P_gamma init_gamma20
#define C_gamma init_gamma20()

#define GEOM_CHECK geom_check20
#endif


#ifdef GEOM11
#define cliff_p 1
#define cliff_q 1

#define P_actionD2 actionD2t11
#define P_actionD4D2 actionD4D2t11
#define P_actionD4D2_r actionD4D2t11_debug
#define C_actionD2 actionD2t11()
#define C_actionD4D2 actionD4D2t11()
#define C_actionD4D2_r actionD4D2t11_debug()

#define P_gamma init_gamma11
#define C_gamma init_gamma11()

#define GEOM_CHECK geom_check11
#endif

#ifdef GEOM01
#define cliff_p 0
#define cliff_q 1

#define P_actionD2 actionD2t01
#define P_actionD4D2 actionD4D2t01
#define C_actionD2 actionD2t01()
#define C_actionD4D2 actionD4D2t01()

#define P_gamma init_gamma01
#define C_gamma init_gamma01()

#define GEOM_CHECK geom_check01
#endif

#ifdef GEOM10
#define cliff_p 1
#define cliff_q 0

#define P_actionD2 actionD2t10
#define P_actionD4D2 actionD4D2t10
#define C_actionD2 actionD2t10()
#define C_actionD4D2 actionD4D2t10()

#define P_gamma init_gamma10
#define C_gamma init_gamma10()

#define GEOM_CHECK geom_check10
#endif

#ifdef GEOM02
#define cliff_p 0
#define cliff_q 2

#define P_actionD2 actionD2t02
#define P_actionD4D2 actionD4D2t02
#define C_actionD2 actionD2t02()
#define C_actionD4D2 actionD4D2t02()

#define P_gamma init_gamma02
#define C_gamma init_gamma02()

#define GEOM_CHECK geom_check02
#endif

#ifdef GEOM03
#define cliff_p 0
#define cliff_q 3

#define P_actionD2 actionD2t03
#define P_actionD4D2 actionD4D2t03
#define C_actionD2 actionD2t03()
#define C_actionD4D2 actionD4D2t03()

#define P_gamma init_gamma03
#define C_gamma init_gamma03()

#define GEOM_CHECK geom_check03
#endif

#ifdef GEOM13
#define cliff_p 1
#define cliff_q 3

#define P_actionD2 actionD2t13
#define P_actionD4D2 actionD4D2t13
#define C_actionD2 actionD2t13()
#define C_actionD4D2 actionD4D2t13()

#define P_gamma init_gamma13
#define C_gamma init_gamma13()

#define GEOM_CHECK geom_check13
#endif

#ifdef GEOM04
#define cliff_p 0
#define cliff_q 4

#define P_actionD2 actionD2t04
#define P_actionD4D2 actionD4D2t04
#define C_actionD2 actionD2t04()
#define C_actionD4D2 actionD4D2t04()

#define P_gamma init_gamma04
#define C_gamma init_gamma04()

#define GEOM_CHECK geom_check04
#endif

#ifdef GEOM40
#define cliff_p 4
#define cliff_q 0

#define P_actionD2 actionD2t40
#define P_actionD4D2 actionD4D2t40
#define C_actionD2 actionD2t40()
#define C_actionD4D2 actionD4D2t40()

#define P_gamma init_gamma40
#define C_gamma init_gamma40()

#define GEOM_CHECK geom_check40
#endif

#ifdef GEOM22
#define cliff_p 2
#define cliff_q 2

#define P_actionD2 actionD2t22
#define P_actionD4D2 actionD4D2t22
#define C_actionD2 actionD2t22()
#define C_actionD4D2 actionD4D2t22()

#define P_gamma init_gamma22
#define C_gamma init_gamma22()

#define GEOM_CHECK geom_check22
#endif

#ifdef GEOM31
#define cliff_p 3
#define cliff_q 1

#define P_actionD2 actionD2t31
#define P_actionD4D2 actionD4D2t31
#define C_actionD2 actionD2t31()
#define C_actionD4D2 actionD4D2t31()

#define P_gamma init_gamma31
#define C_gamma init_gamma31()

#define GEOM_CHECK geom_check31
#endif




#endif
