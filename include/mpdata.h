#ifndef _H_MPDATA
#define _H_MPDATA

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include "./mathfuncs.h"
#include "./vector.h"

#define MPDATA_CYCLIC_X 0x0001U
#define MPDATA_CONTIN_X 0x0002U
#define MPDATA_NILL_X   0x0004U
#define MPDATA_CYCLIC_Y 0x0010U
#define MPDATA_CONTIN_Y 0x0020U
#define MPDATA_NILL_Y   0x0040U
#define MPDATA_CYCLIC_Z 0x0100U
#define MPDATA_CONTIN_Z 0x0200U
#define MPDATA_NILL_Z   0x0400U
#define MPDATA_NONOSCIL 0x1000U

#define MPDATA_OPTS_STD (MPDATA_NILL_X | MPDATA_NILL_Y | MPDATA_NILL_Z | MPDATA_NONOSCIL)

_CCODE_START
typedef struct 
{
    double *data;
    double *G;
    double *u;
    long N;
    long lh;
} arakawa_1d;

typedef struct
{
    double **data;
    double **G;
    double **u[2];
    long N, M;
    long lh;
} arakawa_2d;


#define ARA1D_SZ(GRID) (GRID)->N
#define ARA2D_SZ(GRID) ((GRID)->N * (GRID)->M)
#define ARA1D_HALO(GRID) ((GRID)->data - (GRID)->lh)
#define ARA1D_D(GRID) (GRID)->data

// 1D
arakawa_1d* make_arakawa_1d(vec* data, vec* x, vec* u, vec* G, const double dt);
void memcpy_arakawa_1d(vec* data, arakawa_1d* grid);
vec* get_vec_arakawa(arakawa_1d* grid);
vec* free_arakawa_1d(arakawa_1d* grid, const int get_field);

arakawa_1d* mpadvec1d(arakawa_1d* grid, const long order, const unsigned options);


// 2D
arakawa_2d* make_arakawa_2d(mat* data, vec* x, vec* y, mat* u, mat* v, mat* G, const double dt);
void memcpy_arakawa_2d(mat* data, arakawa_2d* grid);
mat* get_mat_arakawa(arakawa_2d* grid);
mat* free_arakawa_2d(arakawa_2d* grid, const int get_field);

arakawa_2d* mpadvec2d(arakawa_2d* grid, const long order, const unsigned options);

_CCODE_END

#endif