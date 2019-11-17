#include "../include/environment.h"

/*Droplet Advection Physics*/
vec* courant_bin(vec* drdt, vec* dr, const double S, const double dt)
{
    assert(LEN(drdt) == LEN(dr));
    vec* C = allocate_vec(LEN(dr));
    unsigned long i;
    for (i = 0; i < LEN(C); ++i) C->data[i] = S * drdt->data[i] * dt / dr->data[i];
    return C;
}

vec* advec_bin(vec* b, vec* C_drdt)
{
    assert(LEN(b) == LEN(C_drdt));
    vec* fp = zero_vec(LEN(b));
    vec* fm = zero_vec(LEN(b));
    vec* w  = zero_vec(LEN(b));
    double cl, cr;
    double x1, x2;
    double a[3];
    unsigned long i;

    cr = MAX(C_drdt->data[0], 0.0);
    fp->data[0] = MIN(b->data[0], cr * (b->data[0] + (1. - cr) * (b->data[1] - b->data[0]) * 0.5));
    w->data[0] = 1.0;

    for (i = 1; i < LEN(b) - 1; ++i)
    {
        a[0] = (26. * b->data[i] - b->data[i+1] - b->data[i-1]) / 24.;
        a[1] = (b->data[i+1] - b->data[i-1]) / 16.;
        a[2] = (b->data[i+1] - 2. * b->data[i] + b->data[i-1]) /48.;

        cl = -MIN(0., C_drdt->data[i-1]);
        x1 = 1. - 2 * cl;
        x2 = x1 * x1;
        fm->data[i-1] = MAX(0., a[0] * cl - a[1] * (1. - x2) + a[2] * (1. - x1 * x2));

        cr = MAX(0.0, C_drdt->data[i]);
        x1 = 1. - 2. * cr;
        x2 = x1 * x1;
        fp->data[i] = MAX(0., a[0] * cr + a[1] * (1. - x2) + a[2] * (1. - x1 * x2));

        w->data[i] = b->data[i] / MAX(fm->data[i-1] + fp->data[i] + 1e-15, a[0] + 2 * a[2]);
    }

    cl = -MIN(0.0, C_drdt->data[i-1]);
    fm->data[i-1] = MIN(b->data[i], cl * (b->data[i] - (1. - cl) * (b->data[i] - b->data[i-1]) * 0.5));
    w->data[i] = 1.0;

    vec* bn = copy_vec(b);
    for (i = 1; i < b->N - 1; ++i)
        bn->data[i] = bn->data[i] - (fm->data[i-1] + fp->data[i]) * w->data[i] +\
                                     fm->data[i] * w->data[i+1] +\
                                     fp->data[i-1] * w->data[i-1];

    free_vec(fp); free_vec(fm); free_vec(w);

    return bn;
}