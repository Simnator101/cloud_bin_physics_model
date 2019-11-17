#include "../include/mathfuncs.h"

double va_max(unsigned long N, ...)
{
    if (N == 0) return 0.0;
    va_list va;
    unsigned long i;
    double res = -INFINITY;
    double val = 0.0;
    va_start(va, N);

    for (i = 0; i < N; ++i)
    {
        val = va_arg(va, double);
        res = MAX(res, val);
    }
    va_end(va);
    return res;
}

double va_min(unsigned long N, ...)
{
    if (N == 0) return 0.0;
    va_list va;
    unsigned long i;
    double res = INFINITY;
    double val = 0.0;
    va_start(va, N);

    for (i = 0; i < N; ++i)
    {
        val = va_arg(va, double);
        res = MIN(res, val);
    }
    va_end(va);
    return res;
}

double trapz(vec* y, vec* x)
{
    assert(y != NULL && x != NULL);
    assert(y->N == x->N);
    //if (y->N == 1) return 0.0;

    double dx, res = 0.0;
    unsigned long i = 0;

    for (i = 1; i < y->N; ++i)
    {
        dx = x->data[i] - x->data[i-1];
        res += (y->data[i-1] + y->data[i]) * dx / 2.0;
    }
    // Say dx_n = dx_[n-1]
    // res += y->data[y->N - 1] * dx / 2.0;
    return res;
}

vec* gradient_vec(vec* v, vec* x)
{
    assert(v->N == x->N);
    vec* dvdx = zero_vec(v->N);
    unsigned long i;
    double hs, hd;
    
    // Left Wing
    hd = x->data[1] - x->data[0];
    dvdx->data[0] = (v->data[1] - v->data[0]) / hd;

    // Centre
    for (i = 1; i < v->N - 1; ++i)
    {
        hd = x->data[i+1] - x->data[i];
        hs = x->data[i] - x->data[i-1];
        dvdx->data[i]  = hs * hs * v->data[i+1];
        dvdx->data[i] += (hd * hd - hs * hs) * v->data[i];
        dvdx->data[i] -= hd * hd * v->data[i-1];
        dvdx->data[i] /= hs * hd * (hd + hs);
    }

    // Right Wing
    hs = x->data[x->N-1] - x->data[x->N-2];
    dvdx->data[x->N-1] = (v->data[x->N-1] - v->data[x->N-2]) / hd;

    return dvdx;
}

mat* gradient_mat(mat* m, vec* xy, const int axis)
{
    assert(axis == 0 || axis == 1);
    mat* dmdxy = allocate_matNxM(ROWS(m), COLS(m));
    unsigned long i;

    if (axis == 1)
    {
        for (i = 0; i < ROWS(m); ++i)
        {
            vec tmp = {COLS(m), m->data + i * COLS(m)};
            vec* out = gradient_vec(&tmp, xy);
            set_mat_row(dmdxy, out, i);
            free_vec(out);
        }
    }
    else
    {
        for (i = 0; i < COLS(m); ++i)
        {
            vec* tmp = slice_mat_col(m, i);
            vec* out = gradient_vec(tmp, xy);
            set_mat_col(dmdxy, out, i);
            free_vec(out);
            free_vec(tmp);
        }
    }
    

    return dmdxy;
}