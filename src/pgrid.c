#include "../include/pgrid.h"

long long ACTIVE_PGRIDS = 0;

pgrid* allocate_pgrid(unsigned long NZ, unsigned long NY, unsigned NX)
{
    if (NX == 0 || NY == 0 || NZ == 0) return NULL;
    pgrid* grid = malloc(sizeof(pgrid));
    grid->data = malloc(sizeof(double) * NZ * NY * NX);

    grid->NZ = NZ;
    grid->NY = NY;
    grid->NX = NX;
    
    ++ACTIVE_PGRIDS;
    return grid;
}

void free_pgrid(pgrid *grid)
{
    if (grid == NULL) return;
    free(grid->data);
    free(grid);
    --ACTIVE_PGRIDS;
}

pgrid* copy_pgrid(pgrid* p)
{
    if (p == NULL) return NULL;
    pgrid* pn = allocate_pgrid(p->NZ, p->NY, p->NX);
    memcpy(pn->data, p->data, sizeof(double) * LEN_PGRID(p));
    return pn;
}

pgrid* fill_pgrid(unsigned long NZ, unsigned long NY, unsigned long NX, double val)
{
    pgrid* pn = allocate_pgrid(NZ, NY, NX);
    unsigned long i;
    for (i = 0; i < LEN_PGRID(pn); ++i)
        pn->data[i] = val;
    return pn;
}

vec* flatten_pgrid(pgrid* p)
{
    assert(p != NULL);
    vec* v = allocate_vec(LEN_PGRID(p));
    memcpy(v->data, p->data, sizeof(double) * LEN_PGRID(p));
    return v;
}

mat* slice_pgrid_z(pgrid* p, unsigned long zi)
{
    assert(zi < p->NZ);
    mat* m = allocate_matNxM(p->NY, p->NX);
    unsigned long i, j;
    for (i = 0; i < PNY(p); ++i)
        for (j = 0; j < PNX(p); ++j)
            m->data[i * m->M + j] = p->data[zi * p->NY * p->NX + i * p->NX + j];

    return m;
}

mat* slice_pgrid_y(pgrid* p, unsigned long yi)
{
    assert(yi < p->NY);
    mat* m = allocate_matNxM(p->NZ, p->NX);
    unsigned long i, j;
    for (i = 0; i < PNZ(p); ++i)
        for (j = 0; j < PNX(p); ++j)
            m->data[i * m->M + j] = p->data[i * p->NY * p->NX + yi * p->NX + j];

    return m;
}

mat* slice_pgrid_x(pgrid* p, unsigned long xi)
{
    assert(xi < p->NX);
    mat* m = allocate_matNxM(p->NZ, p->NY);
    unsigned long i, j;
    for (i = 0; i < PNZ(p); ++i)
        for (j = 0; j < PNY(p); ++j)
            m->data[i * m->M + j] = p->data[i * p->NY * p->NX + j * p->NX + xi];

    return m;
}

void set_pgrid_z(pgrid* p, mat* m, unsigned long zi)
{
    assert(zi < p->NZ);
    assert(p->NY == m->N && p->NX == m->M);
    unsigned long i;
    for (i = 0; i < LEN_MAT(m); ++i)
        p->data[zi * p->NY * p->NX + i] = m->data[i];
}

void set_pgrid_y(pgrid* p, mat* m, unsigned long yi)
{
    assert(yi < p->NY);
    assert(p->NZ == m->N);
    assert(p->NX == m->M);
    unsigned long i, j, c;
    for (i = 0; i < m->N; ++i) for (j = 0; j < m->M; ++j)
    {
        c = i * p->NY * p->NX + yi * p->NX + j;
        p->data[c] = m->data[i * m->M + j];
    }
}

void set_pgrid_x(pgrid* p, mat* m, unsigned long xi)
{
    assert(xi < p->NX);
    assert(p->NZ == m->N);
    assert(p->NY == m->M);
    unsigned long i, j, c;
    for (i = 0; i < m->N; ++i) for (j = 0; j < m->M; ++j)
    {
        c = i * p->NY * p->NX + j * p->NX + xi;
        p->data[c] = m->data[i * m->M + j];
    }
}

vec* slice_pgrid_zy(pgrid* p, unsigned long zi, unsigned long yi)
{
    assert(zi < p->NZ);
    assert(yi < p->NY);
    vec* v = allocate_vec(p->NX);
    memcpy(v->data, p->data + zi * p->NY * p->NX + yi * p->NX, v->N * sizeof(double));
    return v;
}

vec* slice_pgrid_yx(pgrid* p, unsigned long yi, unsigned long xi)
{
    assert(yi < p->NY);
    assert(xi < p->NX);
    vec* v = allocate_vec(p->NZ);
    unsigned long i;
    for (i = 0; i < PNZ(p); ++i)
        v->data[i] = p->data[i * p->NY * p->NX + yi * p->NX + xi];
    return v;
}

vec* slice_pgrid_zx(pgrid* p, unsigned long zi, unsigned long xi)
{
    assert(zi < p->NZ);
    assert(xi < p->NX);
    vec* v = allocate_vec(p->NY);
    unsigned long i;
    for (i = 0; i < PNY(p); ++i)
        v->data[i] = p->data[zi * p->NY * p->NX + i * p->NX + xi];
    return v;
}

void set_pgrid_zy(pgrid* p, vec* m, unsigned long zi, unsigned long yi)
{
    assert(p->NX == LEN(m));
    unsigned long i, c;
    for (i = 0; i < LEN(m); ++i)
    {
        c = zi * p->NY * p->NX + yi * p->NX + i;
        p->data[c] = m->data[i];
    }
}

void set_pgrid_yx(pgrid* p, vec* m, unsigned long yi, unsigned long xi)
{
    assert(p->NZ == LEN(m));
    unsigned long i, c;
    for (i = 0; i < LEN(m); ++i)
    {
        c = i * p->NY * p->NX + yi * p->NX + xi;
        p->data[c] = m->data[i];
    }
}

void set_pgrid_zx(pgrid* p, vec* m, unsigned long zi, unsigned long xi)
{
    assert(p->NY == LEN(m));
    unsigned long i, c;
    for (i = 0; i < LEN(m); ++i)
    {
        c = zi * p->NY * p->NX + i * p->NX + xi;
        p->data[c] = m->data[i];
    }
}

pgrid* add_pgrid(pgrid* p, unsigned long N, ...)
{
    assert(p != NULL);
    va_list va;
    unsigned long psz = LEN_PGRID(p);
    unsigned long c, i;

    for (c = 0; c < N; ++c)
    {
        pgrid* pv = va_arg(va, pgrid*);
        assert(LEN_PGRID(p) == psz);
        for (i = 0; i < psz; ++i)
            p->data[i] += pv->data[i];
    }
    
    va_start(va, N);
    va_end(va);
    return p;
}


pgrid* sub_pgrid(pgrid* p, unsigned long N, ...)
{
    assert(p != NULL);
    va_list va;
    unsigned long psz = LEN_PGRID(p);
    unsigned long c, i;

    for (c = 0; c < N; ++c)
    {
        pgrid* pv = va_arg(va, pgrid*);
        assert(LEN_PGRID(p) == psz);
        for (i = 0; i < psz; ++i)
            p->data[i] -= pv->data[i];
    }
    
    va_start(va, N);
    va_end(va);
    return p;
}

pgrid* mlt_pgrid(pgrid* p, unsigned long N, ...)
{
    assert(p != NULL);
    va_list va;
    unsigned long psz = LEN_PGRID(p);
    unsigned long c, i;

    for (c = 0; c < N; ++c)
    {
        pgrid* pv = va_arg(va, pgrid*);
        assert(LEN_PGRID(p) == psz);
        for (i = 0; i < psz; ++i)
            p->data[i] *= pv->data[i];
    }
    
    va_start(va, N);
    va_end(va);
    return p;
}

pgrid* div_pgrid(pgrid* p, unsigned long N, ...)
{
    assert(p != NULL);
    va_list va;
    unsigned long psz = LEN_PGRID(p);
    unsigned long c, i;

    for (c = 0; c < N; ++c)
    {
        pgrid* pv = va_arg(va, pgrid*);
        assert(LEN_PGRID(p) == psz);
        for (i = 0; i < psz; ++i)
            p->data[i] /= pv->data[i];
    }
    
    va_start(va, N);
    va_end(va);
    return p;
}

pgrid* scl_pgird(pgrid* p, double scalar)
{
    assert(p != NULL);
    unsigned long i;
    for (i = 0; i < LEN_PGRID(p); ++p)
        p->data[i] *= scalar;
    return p;
}


pgrid* add_scalar_pgrid(pgrid* p, double scalar)
{
    assert(p != NULL);
    unsigned long i;
    for (i = 0; i < LEN_PGRID(p); ++p)
        p->data[i] += scalar;
    return p;
}

pgrid* pow_prid(pgrid* p, double scalar)
{
    assert(p != NULL);
    unsigned long i;
    for (i = 0; i < LEN_PGRID(p); ++p)
        p->data[i] = pow(p->data[i], scalar);
    return p;
}

pgrid* pgrid_map(pgrid* p, double(*fnc)(double))
{
    assert(p != NULL);
    pgrid* pn = allocate_pgrid(p->NZ, p->NY, p->NX);
    unsigned long i; 
    for (i = 0; i < LEN_PGRID(p); ++i) pn->data[i] = (*fnc)(p->data[i]);
    return pn;
}

double pgrid_el_func(pgrid* p, double(*fnc)(double))
{
    assert (p != NULL);
    double res = 0.0;
    unsigned long i;
    for (i = 0; i < LEN_PGRID(p); ++i)
        res += (*fnc)(p->data[i]);
    return res;
}