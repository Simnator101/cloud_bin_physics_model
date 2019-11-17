#include "../include/matrix.h"

long long ACTIVE_MATRICES = 0;

mat* allocate_matNxM(const unsigned long N, const unsigned long M)
{
    assert(N > 0); assert(M > 0);
    mat* r = malloc(sizeof(mat));

    r->data = malloc(sizeof(double) * N * M);
    r->N = N;
    r->M = M;

    ++ACTIVE_MATRICES;
    return r;
}

void free_mat(mat* m)
{
    if (m == NULL) return;
    free(m->data);
    free(m);
    --ACTIVE_MATRICES;
}

mat* copy_mat(mat* m)
{
    mat* res = allocate_matNxM(m->N, m->M);
    memcpy(res->data, m->data, sizeof(double) * m->N * m->M);
    return res;
}

mat* fill_mat(unsigned long N, unsigned M, double val)
{
    mat* m = allocate_matNxM(N, M);
    unsigned i;
    for (i = 0; i < N * M; ++i) m->data[i] = val;
    return m;
}

vec* flatten_mat(mat* m)
{
    assert(m != NULL);
    vec* r = allocate_vec(m->M * m->N);
    memcpy(r->data, m->data, m->N * m->M * sizeof(double));
    return r;
}

mat* reshape_mat(mat* m, unsigned long M, unsigned long N)
{
    assert(M * N == m->N * m->M);
    m->M = M;
    m->N = N;
    return m;
}

mat* reshape_vec(vec* v, unsigned long M, unsigned long N)
{
    assert(M * N == v->N);
    mat* m = malloc(sizeof(mat));

    m->data = malloc(sizeof(double) * v->N);
    memcpy(m->data, v->data, sizeof(double) * M * N);
    m->M = M;
    m->N = N;
    return m;
}

mat* transpose(mat* m)
{
    mat* r = allocate_matNxM(m->M, m->N);
    unsigned long i, j;
    for (i = 0; i < m->N; ++i)
        for (j = 0; j < m->M; ++j)
            r->data[j * m->M + i] = m->data[i * m->M + j];

    free(m->data);
    m->data = r->data;
    m->N = r->N;
    m->M = r->M;
    return m;

}

vec* slice_mat_row(mat* m, unsigned long row)
{
    assert(m != NULL); assert(row < m->N);
    vec* r = allocate_vec(m->M);
    memcpy(r->data, m->data + row * m->M, sizeof(double) * r->N);
    return r;
}

vec* slice_mat_col(mat* m, unsigned long col)
{
    assert(m != NULL); assert(col < m->M);
    vec* r = allocate_vec(m->N);
    unsigned long i;
    for (i = 0; i < m->N; ++i)
        r->data[i] = m->data[i * m->M + col];
    return r;
}

void set_mat_row(mat* m, vec* v, unsigned long row)
{
    assert(m->M == v->N);
    assert(row < m->N);
    memcpy(m->data + row * m->M, v->data, sizeof(double) * v->N);
}

void set_mat_col(mat* m, vec* v, unsigned long col)
{
    assert(m->N == v->N);
    assert(col < m->M);
    unsigned long i; 
    for (i = 0; i < m->N; ++i)
        m->data[i * m->M + col] = v->data[i];
}

void set_mat_ind(mat* m, unsigned long row, unsigned long col, const double val)
{
    assert(row < m->N);
    assert(col < m->M);
    m->data[row * m->M + col] = val;
}

vec __vec_wrapper_mat(mat* m)
{
    vec v = {m->M * m->N, m->data};
    return v;
}

mat* add_mat(mat* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long matsz = LEN_MAT(to);
    unsigned long c, i; 

    va_start(va, N);
    for (c = 0; c < N; ++c)
    {
        mat* m = va_arg(va, mat*);
        assert(LEN_MAT(m) == matsz);
        for (i = 0; i < matsz; ++i)
            to->data[i] += m->data[i];
    }

    va_end(va);
    return to;
}

mat* sub_mat(mat* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long matsz = LEN_MAT(to);
    unsigned long c, i; 

    va_start(va, N);
    for (c = 0; c < N; ++c)
    {
        mat* m = va_arg(va, mat*);
        assert(LEN_MAT(m) == matsz);
        for (i = 0; i < matsz; ++i)
            to->data[i] -= m->data[i];
    }

    va_end(va);
    return to;
}

mat* mlt_mat(mat* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long matsz = LEN_MAT(to);
    unsigned long c, i; 

    va_start(va, N);
    for (c = 0; c < N; ++c)
    {
        mat* m = va_arg(va, mat*);
        assert(LEN_MAT(m) == matsz);
        for (i = 0; i < matsz; ++i)
            to->data[i] *= m->data[i];
    }

    va_end(va);
    return to;
}

mat* div_mat(mat* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long matsz = LEN_MAT(to);
    unsigned long c, i; 

    va_start(va, N);
    for (c = 0; c < N; ++c)
    {
        mat* m = va_arg(va, mat*);
        assert(LEN_MAT(m) == matsz);
        for (i = 0; i < matsz; ++i)
            to->data[i] /= m->data[i];
    }

    va_end(va);
    return to;
}

mat* scl_mat(mat* to, double scaler)
{
    vec __tmpv = __vec_wrapper_mat(to);
    scl_vec(&__tmpv, scaler);
    return to;
}


mat* add_scalar_mat(mat* to, double scalar)
{
    vec __tmpv = __vec_wrapper_mat(to);
    add_scalar_vec(&__tmpv, scalar);
    return to;
}

mat* pow_mat(mat* to, double scalar)
{
    vec __tmpv = __vec_wrapper_mat(to);
    pow_vec(&__tmpv, scalar);
    return to;
}