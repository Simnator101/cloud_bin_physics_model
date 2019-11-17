#ifndef _H_MATRIX
#define _H_MATRIX

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

#include "./vector.h"

_CCODE_START
typedef struct __matrix
{
    unsigned long N, M;
    double* data;
} mat;

extern long long ACTIVE_MATRICES;

#define LEN_MAT(D) ((D)->M * (D)->N)
#define ROWS(D) ((D)->N)
#define COLS(D) ((D)->M)


mat* allocate_matNxM(const unsigned long N, const unsigned long M);
#define allocate_matNxN(N) allocate_matNxM(N, N)
void free_mat(mat* m);

mat* copy_mat(mat* m);
mat* fill_mat(unsigned long N, unsigned M, double val);
#define zero_mat(N,M) fill_mat((N), (M), 0.0)
#define ones_mat(N,M) fill_mat((N), (M), 1.0)

vec* flatten_mat(mat* m);
mat* transpose(mat* m);
mat* reshape_mat(mat* m, unsigned long M, unsigned long N);
mat* reshape_vec(vec* v, unsigned long M, unsigned long N);

vec* slice_mat_row(mat* m, unsigned long row);
vec* slice_mat_col(mat* m, unsigned long col);
void set_mat_row(mat* m, vec* v, unsigned long row);
void set_mat_col(mat* m, vec* v, unsigned long col);
void set_mat_ind(mat* m, unsigned long row, unsigned long col, const double val);

mat* add_mat(mat* to, const unsigned long N, ...);
mat* sub_mat(mat* to, const unsigned long N, ...);
mat* mlt_mat(mat* to, const unsigned long N, ...);
mat* div_mat(mat* to, const unsigned long N, ...);

mat* scl_mat(mat* to, double scaler);
mat* add_scalar_mat(mat* to, double scalar);
mat* pow_mat(mat* to, double scalar);

#define max_mat(v) __max_in_buffer((v)->data, LEN_MAT(v))
#define min_mat(v) __min_in_buffer((v)->data, LEN_MAT(v))

void fprint_mat(FILE* pf, mat* m);

_CCODE_END


#endif