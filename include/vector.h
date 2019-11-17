#ifndef _H_VECTOR
#define _H_VECTOR

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <memory.h>
#include <assert.h>

#include <math.h>

#include "constants.h"

_CCODE_START

double __min_in_buffer(double* data, unsigned long els);
double __max_in_buffer(double* data, unsigned long els);


typedef struct __vector
{
    unsigned long N;
    double* data;
} vec;

extern long long ACTIVE_VECTORS;

#define LEN(V) (V->N)


vec* allocate_vec(const unsigned long N);
void free_vec(vec* v);

vec* copy_vec(vec* v);
vec* fill_vec(const unsigned long N, double val);
#define zero_vec(N) fill_vec(N, 0.0)
#define one_vec(N) fill_vec(N, 1.0)
vec* range(double start, double end, double step_sz);
vec* linspace(double start, double end, unsigned long N);

vec* push_back(vec* v, double val);
vec* push_front(vec* v, double val);
vec* slice_vec(vec* iv, const unsigned long left, const unsigned long right);
void set_vec(vec* iv, vec* v);
vec* pop_back(vec* v);
vec* pop_front(vec* v);

vec* add_vec(vec* to, const unsigned long N, ...);
vec* sub_vec(vec* to, const unsigned long N, ...);
vec* mlt_vec(vec* to, const unsigned long N, ...);
vec* div_vec(vec* to, const unsigned long N, ...);

vec* scl_vec(vec* to, double scaler);
vec* add_scalar_vec(vec* to, double scalar);
vec* pow_vec(vec* to, double scalar);
vec* clip_vec(vec* to, double lower, double upper);

double sum_vec(vec* v);
double prod_vec(vec* v);
double dot_vec(vec* a, vec* b);
#define max_vec(v) __max_in_buffer((v)->data, LEN(v))
#define min_vec(v) __min_in_buffer((v)->data, LEN(v))


double interp(double x, vec* x0, vec* y0, interp_type mode);
vec* interp_vec(vec* x, vec* x0, vec* y0, interp_type mode);

void fprint_vec(FILE* pf, vec* v);

_CCODE_END

#endif