#ifndef _H_MATHFUNCS
#define _H_MATHFUNCS

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include <stdarg.h>
#include <math.h>
#include "./vector.h"
#include "./matrix.h"

#ifndef MAX
#define MAX(A,B) ((A) >= (B) ? (A) : (B))
#endif
#ifndef MIN
#define MIN(A,B) ((A) <= (B) ? (A) : (B))
#endif

#define BOUNDARY_NILL 0x01
#define BOUNDARY_CYCLIC 0x02
#define BOUNDARY_CONTINUOUS 0x04
#define BOUNDARY_KEEP_OUT_FLUXES 0x08

_CCODE_START
/*va arg reduction*/
double va_max(unsigned long N, ...);
double va_min(unsigned long N, ...);

/*Integration*/
double trapz(vec* y, vec* x);

/*Gradients*/
vec* gradient_vec(vec* v, vec* x);
mat* gradient_mat(mat* m, vec* xy, const int axis);

/*Definite Positive Advection Routines by Andreas Bott*/
//vec* courant_bott(const vec* v, const vec* x, const double dt);
//vec* adv2p(vec* y, vec* c, vec* rho, int boundary);

/*MPDATA Scheme*/
//mat* mpadvec(mat* phi, mat* cu, mat* cv, const mat* G, int order, int boundary);

_CCODE_END
#endif