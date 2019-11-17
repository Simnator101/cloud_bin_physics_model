#ifndef _H_PGRID
#define _H_PGRID

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

#include "./vector.h"
#include "./matrix.h"

_CCODE_START
typedef struct __pgrid
{
    unsigned long NZ, NY, NX;
    double *data;
} pgrid;

extern long long ACTIVE_PGRIDS;

#define LEN_PGRID(G) ((G)->NZ * (G)->NY * (G)->NX)
#define PNX(G) (G)->NX
#define PNY(G) (G)->NY
#define PNZ(G) (G)->NZ

pgrid* allocate_pgrid(unsigned long NZ, unsigned long NY, unsigned NX);
void free_pgrid(pgrid *grid);

pgrid* copy_pgrid(pgrid* p);
pgrid* fill_pgrid(unsigned long NZ, unsigned long NY, unsigned long NX, double val);
#define zero_pgrid(NZ,NY,NX) fill_pgrid((NZ), (NY), (NX), 0.0)
#define ones_pgrid(NZ,NY,NX) fill_pgrid((NZ), (NY), (NX), 1.0)

vec* flatten_pgrid(pgrid* p);

mat* slice_pgrid_z(pgrid* p, unsigned long zi);
mat* slice_pgrid_y(pgrid* p, unsigned long yi);
mat* slice_pgrid_x(pgrid* p, unsigned long xi);

void set_pgrid_z(pgrid* p, mat* m, unsigned long zi);
void set_pgrid_y(pgrid* p, mat* m, unsigned long yi);
void set_pgrid_x(pgrid* p, mat* m, unsigned long xi);

vec* slice_pgrid_zy(pgrid* p, unsigned long zi, unsigned long yi);
vec* slice_pgrid_yx(pgrid* p, unsigned long yi, unsigned long xi);
vec* slice_pgrid_zx(pgrid* p, unsigned long zi, unsigned long xi);

void set_pgrid_zy(pgrid* p, vec* m, unsigned long zi, unsigned long yi);
void set_pgrid_yx(pgrid* p, vec* m, unsigned long yi, unsigned long xi);
void set_pgrid_zx(pgrid* p, vec* m, unsigned long zi, unsigned long xi);


pgrid* add_pgrid(pgrid* p, unsigned long N, ...);
pgrid* sub_pgrid(pgrid* p, unsigned long N, ...);
pgrid* mlt_pgrid(pgrid* p, unsigned long N, ...);
pgrid* div_pgrid(pgrid* p, unsigned long N, ...);

pgrid* scl_pgird(pgrid* p, double scalar);
pgrid* add_scalar_pgrid(pgrid* p, double scalar);
pgrid* pow_prid(pgrid* p, double scalar);

#define max_pgrid(v) __max_in_buffer((v)->data, LEN_PGRID(v))
#define min_pgrid(v) __min_in_buffer((v)->data, LEN_PGRID(v))

pgrid* pgrid_map(pgrid* p, double(*fnc)(double));
double pgrid_el_func(pgrid* p, double(*fnc)(double));

_CCODE_END

#endif