#include "../include/vector.h"

long long ACTIVE_VECTORS = 0;

double __min_in_buffer(double* data, unsigned long els)
{
    double min = 1e100;
    unsigned long i;
    for (i = 0; i < els; ++i)
        min = (data[i] < min) ? data[i] : min;
    return min;
}

double __max_in_buffer(double* data, unsigned long els)
{
    double max = -1e100;
    unsigned long i;
    for (i = 0; i < els; ++i)
        max = (data[i] > max) ? data[i] : max;
    return max;
}

vec* allocate_vec(const unsigned long N)
{
    vec* v = (vec*)malloc(sizeof(vec));
    v->data = (N > 0) ? malloc(N * sizeof(double)) : NULL;
    v->N = N;
    ++ACTIVE_VECTORS;
    return v;
}

void free_vec(vec* v)
{
    if (v == NULL) return;
    free(v->data);
    free(v);
    ACTIVE_VECTORS = (ACTIVE_VECTORS > 0) ? ACTIVE_VECTORS - 1 : 0;
}


vec* copy_vec(vec* v)
{
    if (v == NULL) return NULL;
    vec* cv = allocate_vec(v->N);
    memcpy(cv->data, v->data, v->N * sizeof(double));
    return cv;
}

vec* fill_vec(const unsigned long N, double val)
{
    unsigned long i;
    vec* v = allocate_vec(N);
    for (i = 0; i < N; ++i) v->data[i] = val;
    return v;
}

vec* range(double start, double end, double step_sz)
{
    unsigned long c = 0, i;
    vec* v = NULL;
    assert(start != end);
    if (end < start) assert(step_sz < 0.0);
    if (end > start) assert(step_sz > 0.0);

    while(start + step_sz * (double)c < end) ++c;
    v = allocate_vec(c);
    for (i = 0; i < c; ++i)
        v->data[i] = start + (double)i * step_sz;

    return v;
}

vec* linspace(double start, double end, unsigned long N)
{
    if (N <= 1) return NULL;
    assert(end != start);

    unsigned long i;
    double dydx = (end - start) / ((double)N - 1);
    vec* v = allocate_vec(N);
    for (i = 0; i < N; ++i)
        v->data[i] = start + (double)i * dydx;
    return v;
}


vec* push_back(vec* v, double val)
{
    if (v == NULL) return fill_vec(1, val);
    v->data = realloc(v->data, (++v->N) * sizeof(double));
    v->data[v->N-1] = val;
    return v;
}

vec* push_front(vec* v, double val)
{
    if (v == NULL) return fill_vec(1, val);
    double* dnew = malloc((++v->N) * sizeof(double));
    memcpy(dnew + 1, v->data, sizeof(double) * (v->N - 1));
    dnew[0] = val;
    free(v->data);
    v->data = dnew;
    return v;
}

vec* slice_vec(vec* iv, const unsigned long left, const unsigned long right)
{
    assert(iv != NULL);
    assert(left < right);
    assert(right <= iv->N);

    vec* v = allocate_vec(right - left);
    memcpy(v->data, iv->data + left, (right - left) * sizeof(double));
    return v;
}

void set_vec(vec* iv, vec* v)
{
    assert(iv->N == v->N);
    memcpy(iv->data, v->data, sizeof(double) * v->N);
}

vec* pop_back(vec* v)
{
    if (v->N == 0) return v;
    v->data = realloc(v->data, --v->N * sizeof(double));
    return v;
}

vec* pop_front(vec* v)
{
    if (v->N == 0) return v;
    double* ndata = malloc(--v->N * sizeof(double));
    memcpy(ndata, v->data + 1, sizeof(double) * v->N);
    free(v->data);
    v->data = ndata;
    return v;
}

vec* add_vec(vec* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long i,j;  

    va_start(va, N);
    for (i = 0; i < N; ++i)
    {
        vec* el = va_arg(va, vec*);
        assert(el->N == to->N);
        for (j = 0; j < el->N; ++j)
            to->data[j] += el->data[j];
    }  
    va_end(va);

    return to;
}

vec* sub_vec(vec* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long i,j;  

    va_start(va, N);
    for (i = 0; i < N; ++i)
    {
        vec* el = va_arg(va, vec*);
        assert(el->N == to->N);
        for (j = 0; j < el->N; ++j)
            to->data[j] -= el->data[j];
    }  
    va_end(va);

    return to;
}

vec* mlt_vec(vec* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long i,j;  

    va_start(va, N);
    for (i = 0; i < N; ++i)
    {
        vec* el = va_arg(va, vec*);
        assert(el->N == to->N);
        for (j = 0; j < el->N; ++j)
            to->data[j] *= el->data[j];
    }  
    va_end(va);

    return to;
}

vec* div_vec(vec* to, const unsigned long N, ...)
{
    assert(to != NULL);
    va_list va;
    unsigned long i,j;  

    va_start(va, N);
    for (i = 0; i < N; ++i)
    {
        vec* el = va_arg(va, vec*);
        assert(el->N == to->N);
        for (j = 0; j < el->N; ++j)
            to->data[j] /= el->data[j];
    }  
    va_end(va);

    return to;
}

vec* scl_vec(vec* to, double scaler)
{
    assert(to != NULL);
    unsigned long i;
    for (i = 0; i < to->N; ++i)
        to->data[i] *= scaler;
    return to;
}

vec* add_scalar_vec(vec* to, double scalar)
{
    assert(to != NULL);
    unsigned long i;
    for (i = 0; i < to->N; ++i)
        to->data[i] += scalar;
    return to;
}

vec* pow_vec(vec* to, double scalar)
{
    assert(to != NULL);
    unsigned long i;
    for (i = 0; i < to->N; ++i)
        to->data[i] = pow(to->data[i], scalar);
    return to;
}

vec* clip_vec(vec* to, double lower, double upper)
{
    assert(lower < upper);
    unsigned long i;
    for (i = 0; i < LEN(to); ++i)
    {
        double lv = (to->data[i] > lower) ? to->data[i] : lower;
        double rv = (lv < upper) ? lv : upper;
        to->data[i] = rv;
    }
    return to;
}

double sum_vec(vec* v)
{
    if (v == NULL) return 0.0;
    double sm = 0.0;
    unsigned long i;
    for (i = 0; i < v->N; ++i)
        sm += v->data[i];
    return sm;
}

double prod_vec(vec* v)
{
    if (v == NULL || v->N == 0) return 0.0;
    double prod = v->data[0];
    unsigned long i;
    for (i = 1; i < v->N; ++i)
        prod *= v->data[i];
    return prod;
}

double dot_vec(vec* a, vec* b)
{
    if (a == NULL || b == NULL) return NAN;
    assert(a->N == b->N);
    double dot = 0.0;
    unsigned long i;
    for (i = 0; i < a->N; ++i)
        dot += a->data[i] * b->data[i];
    return dot;
}

double interp(double x, vec* x0, vec* y0, interp_type mode)
{
    assert(x0 != NULL && y0 != NULL);
    assert(x0->N == y0->N);

    if (x >= x0->data[x0->N-1]) return y0->data[y0->N-1];
    if (x <= x0->data[0]) return y0->data[0];

    unsigned li = 0, ri = x0->N - 1;
    double f = 0.0;
    while (x0->data[li + 1] < x) ++li;
    while (ri > 0 && x0->data[ri - 1] > x) --ri;
    f = (x0->data[ri] - x) / (x - x0->data[li] + 1e-40);
    f = 1. / (f + 1.);

    if (mode == LINEAR)
       return f * y0->data[ri] + (1. - f) * y0->data[li];
    else if (mode == LOG)
       return pow(y0->data[ri], f) * pow(y0->data[li], 1. - f);
    return NAN;
}

vec* interp_vec(vec* x, vec* x0, vec* y0, interp_type mode)
{
    assert(x != NULL);
    assert(x0 != NULL && y0 != NULL);
    assert(x0->N == y0->N);

    vec* v = allocate_vec(x->N);
    unsigned long i, ri, li;
    double f;

    for (i = 0; i < x->N; ++i)
    {
        // First two sections for out of bounds
        if (x->data[i] >= x0->data[x0->N-1])
        {
            v->data[i] = y0->data[y0->N-1];
        }
        else if (x->data[i] <= x0->data[0])
        {
            v->data[i] = y0->data[0];
        }
        else
        {
            ri = x0->N-1; li = 0;
            while (x0->data[li + 1] < x->data[i])
                ++li;
            while (ri > 0 && x0->data[ri - 1] > x->data[i])
                --ri;

            f = (x0->data[ri] - x->data[i]) / (x->data[i] - x0->data[li] + 1e-40);
            f = 1. / (f + 1.);

            if (mode == LINEAR)
                v->data[i] = f * y0->data[ri] + (1. - f) * y0->data[li];
            else if (mode == LOG)
                v->data[i] = pow(y0->data[ri], f) * pow(y0->data[li], 1. - f);
            
        }
        
    }
    

    return v;
}

void fprint_vec(FILE* pf, vec* v)
{
    assert(v != NULL);
    assert(pf != NULL);
    unsigned long i; 

    fprintf(pf, "vec([");
    for (i = 0; i < v->N; ++i)
    {
        double val = v->data[i];
        if ((fabs(val) > 1e4 || fabs(val) < 1e-4) && (val != 0.0))
        {
            fprintf(pf, "%4E,", val);
        }
        else
        {
            fprintf(pf, "%4f,", val);
        }

        if (i != 0 && (i % 5) == 0 && (i + 1) < v->N)
            fprintf(pf, "\r\n     ");
        else
            fprintf(pf, " ");
        
    }
    fprintf(pf, "\b\b])\r\n");

}