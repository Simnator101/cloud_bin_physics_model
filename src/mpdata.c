#include "../include/mpdata.h"

#define ARA1D_TOTS(GRID) ((GRID)->N + 2 * (GRID)->lh)
#define ARA2D_TOTS(GRID) (((GRID)->N + 2 * (GRID)->lh) * ((GRID)->M + 2 * (GRID)->lh))
#define CPOS(C) ((C + fabs(C)) * 0.5)
#define CNEG(C) ((C - fabs(C)) * 0.5)
#define EPS 1e-20

arakawa_1d* make_arakawa_1d(vec* data, vec* x, vec* u, vec* G, const double dt)
{
    assert(data->N == x->N); assert(x->N == u->N);
    arakawa_1d* grid = malloc(sizeof(arakawa_1d));
    grid->N = x->N;
    grid->lh = 2;
    unsigned long size = x->N + 2 * grid->lh;
    unsigned long i;

    grid->data = (double*)(calloc(size, sizeof(double))) + grid->lh;
    grid->G = (double*)(calloc(size, sizeof(double))) + grid->lh;
    grid->u = (double*)(calloc(size + 1, sizeof(double))) + grid->lh;

    memcpy(grid->data, data->data, sizeof(double) * ARA1D_SZ(grid));
    if (G == NULL)
        for (i = 0; i < x->N; ++i) grid->G[i] = 1.0;
    else
        memcpy(grid->G, G->data, sizeof(double) * G->N);

    
    double dx = x->data[1] - x->data[0];
    grid->u[0] = u->data[0] * dt / dx;
    for (i = 0; i < grid->N - 1; ++i)
    {
        dx = x->data[i+1] - x->data[i];
        grid->u[i+1] = 0.5 * (u->data[i] + u->data[i+1]) * dt / dx;
    }
    dx = x->data[i] - x->data[i-1];
    grid->u[i+1] = u->data[i] * dt / dx; 

    return grid;
}

void memcpy_arakawa_1d(vec* data, arakawa_1d* grid)
{
    assert(data->N == grid->N);
    memcpy(data->data, grid->data, sizeof(double) * grid->N);
}

vec* get_vec_arakawa(arakawa_1d* grid)
{
    assert(grid != NULL);
    vec* ndata = allocate_vec(grid->N);
    memcpy(ndata->data, grid->data, grid->N * sizeof(double));
    return ndata;
}


vec* free_arakawa_1d(arakawa_1d* grid, const int get_field)
{
    vec* res = (get_field) ? get_vec_arakawa(grid) : NULL;
    free(grid->data - grid->lh);
    free(grid->u - grid->lh);
    free(grid->G - grid->lh);
    free(grid);
    return res;
}

#define CPOS(C) ((C + fabs(C)) * 0.5)
#define CNEG(C) ((C - fabs(C)) * 0.5)


static double flux(double psi_l, double psi_r, double C)
{
    return (psi_l * CPOS(C) + psi_r * CNEG(C));
}

static double mp_A(double psi_p, double psi)
{
    return (psi_p - psi) / (psi_p + psi + EPS);
}

static double mp_B2d(double** d, double** cv, double** ov, int i, int j, int axis)
{
    double cav = 0.25 * (ov[i+1][j+1] + ov[i][j+1] + ov[i+1][j] + ov[i][j]);
    if (axis == 0) // Y is the main axis
    {
        double B[2] = 
        {
            d[i+1][j+1] + d[i][j+1],
            d[i+1][j-1] + d[i][j-1]
        };

        return 0.5 * cav * cv[i][j] * (B[0] - B[1]) / (B[0] + B[1] + EPS); 
    }
    double B[2] = 
    {
        d[i+1][j+1] + d[i][j+1],
        d[i+1][j-1] + d[i][j-1]
    };
    return 0.5 * cav * cv[i][j] * (B[0] - B[1]) / (B[0] + B[1] + EPS); 
}

static int fill_halo_x(double* d, double* u, double* G, unsigned long NX, const int options)
{
    if (options & MPDATA_NILL_X)
    {
        d[-2] = 0.0;    d[NX+1] = 0.0;
        d[-1] = 0.0;    d[NX] = 0.0;

        u[-2] = u[0];   u[NX+2] = u[NX];
        u[-1] = u[0];   u[NX+1] = u[NX];

        G[-2] = 1.0;    G[NX+1] = 1.0;
        G[-1] = 1.0;    G[NX] = 1.0;
        return 1;
    }
    else if (options & MPDATA_CONTIN_X)
    {
        d[-2] = d[0];    d[NX+1] = d[NX-1];
        d[-1] = d[0];    d[NX] = d[NX-1];

        u[-2] = u[0];   u[NX+2] = u[NX];
        u[-1] = u[0];   u[NX+1] = u[NX];

        G[-2] = G[0];    G[NX+1] = G[NX-1];
        G[-1] = G[0];    G[NX] = G[NX-1];
        return 1;
    }
    else if (options & MPDATA_CYCLIC_X)
    {
        d[-2] = d[NX-2];    d[NX+1] = d[1];
        d[-1] = d[NX-1];    d[NX] = d[0];

        u[-2] = u[NX-1];      u[NX+2] = u[1];
        u[-1] = u[NX];      u[NX+1] = u[0];

        G[-2] = G[NX-2];    G[NX+1] = G[1];
        G[-1] = G[NX-1];    G[NX] = G[0];
        return 1;
    }
    
    return 0;
}

static int fill_halo_y(double** d, double** v, double** G, unsigned long NY, unsigned long NX, const int options)
{
    memset(d[-2] - 2, 0, sizeof(double) * (NX + 4));
    memset(d[-1] - 2, 0, sizeof(double) * (NX + 4));
    memset(d[NY] - 2, 0, sizeof(double) * (NX + 4));
    memset(d[NY+1] - 2, 0, sizeof(double) * (NX + 4));

    memset(v[-2] - 2, 0, sizeof(double) * (NX + 4));
    memset(v[-2] - 2, 0, sizeof(double) * (NX + 4));
    memset(v[NY+1] - 2, 0, sizeof(double) * (NX + 4));
    memset(v[NY+2] - 2, 0, sizeof(double) * (NX + 4));

    if (options & MPDATA_NILL_Y)
    {
        long i;
        for (i = 0; i < NX; ++i)
        {
            //d[-2][i] = 0.0; d[NY  ][i] = 0.0;
            //d[-1][i] = 0.0; d[NY+1][i] = 0.0;

            v[-2][i] = v[0][i];
            v[-1][i] = v[0][i];

            v[NY+2][i] = v[NY][i];
            v[NY+1][i] = v[NY][i];


            G[-2][i] = 1.0; G[NY    ][i] = 1.0;
            G[-1][i] = 1.0; G[NY+1  ][i] = 1.0;
        }

        

        return 1;
    }
    else if (options & MPDATA_CONTIN_Y)
    {
        unsigned long i;
        for (i = 0; i < NX; ++i)
        {
            d[-2][i] = d[0][i]; d[NY  ][i] = d[NY-1][i];
            d[-1][i] = d[0][i]; d[NY+1][i] = d[NY-1][i];

            v[-2][i] = v[0][i];
            v[-1][i] = v[0][i];

            v[NY+2][i] = v[NY][i];
            v[NY+1][i] = v[NY][i];

            G[-2][i] = G[0][i]; G[NY    ][i] = G[NY-1][i];
            G[-1][i] = G[0][i]; G[NY+1  ][i] = G[NY-1][i];
        }
        return 1;
    }
    else if (options & MPDATA_CYCLIC_Y)
    {
        unsigned long i;
        for (i = 0; i < NX; ++i)
        {
            d[-2][i] = d[NY-2][i]; d[NY  ][i] = d[0][i];
            d[-1][i] = d[NY-1][i]; d[NY+1][i] = d[1][i];

            v[-2][i] = v[NY-1][i];
            v[-1][i] = v[NY][i];

            v[NY+2][i] = v[1][i];
            v[NY+1][i] = v[0][i];

            G[-2][i] = G[NY-2][i]; G[NY    ][i] = G[0][i];
            G[-1][i] = G[NY-1][i]; G[NY+1  ][i] = G[1][i];
        }
        return 1;
    }
    return 0;
}

static double beta_up1d(double* psi, double* C, double* psi_max, const double *G, int i)
{
    double nom = G[i] * (va_max(4, psi_max[i], psi[i-1], psi[i], psi[i+1]) - psi[i]);
    double den = CPOS(C[i]) * psi[i-1] - CNEG(C[i+1]) * psi[i+1];
    return nom / (den + EPS);
}

static double beta_up2d(double** psi, double** C[2], double** psi_max, double **G, int i, int j)
{
    double nom = G[i][j] * (va_max(6, psi_max[i][j],
                                      psi[i-1][j],
                         psi[i][j-1], psi[i][j], psi[i][j+1],
                                      psi[i+1][j]) - psi[i][j]);


    double den =  (CPOS(C[1][i][j]) * psi[i-1][j] - CNEG(C[1][i+1][j]) * psi[i+1][j]) + \
                  (CPOS(C[0][i][j]) * psi[i][j-1] - CNEG(C[0][i][j+1]) * psi[i][j+1]);
    return nom / (den + EPS);
}

static double beta_dn1d(double* psi, double* C, double* psi_min, const double *G, int i)
{
    double nom = G[i] * (psi[i] - va_min(4, psi_min[i], psi[i-1], psi[i], psi[i+1]));
    double den = CPOS(C[i+1]) * psi[i] - CNEG(C[i]) * psi[i];
    return nom / (den + EPS);
}

static double beta_dn2d(double** psi, double** C[2], double** psi_min, double **G, int i, int j)
{
    double nom = G[i][j] * (psi[i][j] - va_min(6, psi_min[i][j],
                                      psi[i-1][j],
                         psi[i][j-1], psi[i][j], psi[i][j+1],
                                      psi[i+1][j]));


    double den = (CPOS(C[1][i+1][j]) * psi[i][j] - CNEG(C[1][i][j]) * psi[i][j]) + \
                 (CPOS(C[0][i][j+1]) * psi[i][j] - CNEG(C[0][i][j]) * psi[i][j]);
    return nom / (den + EPS);
}

arakawa_1d* mpadvec1d(arakawa_1d* grid, const long order, const unsigned options)
{
    // Parse Options
    assert(fill_halo_x(grid->data, grid->u, grid->G, grid->N, options));
    assert(order >= 1);
     
    unsigned FL = options & MPDATA_NONOSCIL;
    long o, i;

    double* psip = grid->data;

    
    double flow;

    // Create temp fields
    //double* C0 = grid->u;
    double* C = (double*)(calloc(ARA1D_TOTS(grid) + 1, sizeof(double))) + grid->lh;
    const double* G = grid->G;
    double* psin = (double*)(calloc(ARA1D_TOTS(grid), sizeof(double))) + grid->lh;
    double* Cdif = (double*)(calloc(ARA1D_TOTS(grid) + 1, sizeof(double))) + grid->lh;

    memcpy(C - grid->lh, grid->u - grid->lh, (ARA1D_TOTS(grid) + 1) * sizeof(double));

    // FL fields
    double* psi_max = NULL;
    double* psi_min = NULL;
    if (FL)
    {
        psi_max = (double*)(calloc(ARA1D_TOTS(grid), sizeof(double))) + grid->lh;
        psi_min = (double*)(calloc(ARA1D_TOTS(grid), sizeof(double))) + grid->lh;
        psi_max[-1] = va_max(3, psip[-2], psip[-1], psip[0]);
        psi_min[-1] = va_min(3, psip[-2], psip[-1], psip[0]);
        for (i = 0; i < grid->N; ++i)
        {
            psi_max[i] = va_max(3, psip[i-1], psip[i], psip[i+1]);
            psi_min[i] = va_min(3, psip[i-1], psip[i], psip[i+1]);
        }
        psi_max[grid->N] = va_max(3, psip[grid->N-1], psip[grid->N], psip[grid->N+1]);
        psi_min[grid->N] = va_min(3, psip[grid->N-1], psip[grid->N], psip[grid->N+1]);
    }

    for (o = 0; o < order; ++o)
    {
        

        // Apply Flux Scheme
        for (i = 0; i < (long)grid->N; ++i)
        {
            flow = flux(psip[i], psip[i+1], C[i+1]) - flux(psip[i-1], psip[i], C[i]);
            psin[i] = psip[i] - flow / G[i];
        }

        // Calculate antidiffusive velocities
        fill_halo_x(psin, C, grid->G, grid->N, options);
        for (i = -1; i < (long)grid->N; ++i)
        {
            Cdif[i+1] = fabs(C[i+1]) * (1.0 - fabs(C[i+1]));
            Cdif[i+1] = Cdif[i+1] * mp_A(psin[i+1], psin[i]);
        }

        // Apply Non Oscillatory Options
        if (!FL) goto skip_oscil;

        double* Cm = (double*)(calloc(ARA1D_TOTS(grid) + 1, sizeof(double))) + grid->lh;
        for (i = -1; i < (long)grid->N; ++i)
        {
            double bscale = 1.0;
            if (Cdif[i+1] >= 0.0)
                if (psin[i] > 0.0)
                    bscale = va_min(3, 1.0, beta_dn1d(psin, Cdif, psi_min, G, i), beta_up1d(psin, Cdif, psi_max, G, i+1));
                else
                    bscale = va_min(3, 1.0, beta_up1d(psin, Cdif, psi_max, G, i), beta_dn1d(psin, Cdif, psi_min, G, i+1));             
            else
                if (psin[i+1] > 0.0)
                    bscale = va_min(3, 1.0, beta_up1d(psin, Cdif, psi_max, G, i), beta_dn1d(psin, Cdif, psi_min, G, i+1));
                else
                    bscale = va_min(3, 1.0, beta_dn1d(psin, Cdif, psi_min, G, i), beta_up1d(psin, Cdif, psi_max, G, i+1));
            Cm[i+1] = Cdif[i+1] * bscale;
        }
        memcpy(Cdif - grid->lh, Cm - grid->lh, (ARA1D_TOTS(grid) + 1) * sizeof(double));
        free(Cm - grid->lh);

        skip_oscil:
        fill_halo_x(psin, Cdif, grid->G, grid->N, options);
        memcpy(psip - grid->lh, psin - grid->lh, sizeof(double) * ARA1D_TOTS(grid));
        memcpy(C - grid->lh, Cdif - grid->lh, (ARA1D_TOTS(grid) + 1) * sizeof(double));
        //C = Cdif;
    }
    free(psin - grid->lh);
    free(Cdif - grid->lh);
    free(C - grid->lh);
    if (psi_max != NULL) free(psi_max - grid->lh);
    if (psi_min != NULL) free(psi_min - grid->lh);
    return grid;
}

arakawa_2d* make_arakawa_2d(mat* data, vec* x, vec* y, mat* u, mat* v, mat* G, const double dt)
{
    assert((ROWS(data) == LEN(y)) && (COLS(data) == LEN(x)));
    assert(ROWS(data) == ROWS(u));
    assert(ROWS(u) == ROWS(v));
    assert(COLS(data) == COLS(u));
    assert(COLS(u) == COLS(v));

    arakawa_2d* grid = malloc(sizeof(arakawa_2d));
    grid->N = LEN(y); grid->M = LEN(x);
    grid->lh = 2;

    unsigned long i,j;
    // Allocate data matrix with halo
    grid->data = calloc(2 * grid->lh + ROWS(data), sizeof(double*));
    for (i = 0; i < ROWS(data) + 2 * grid->lh; ++i) 
    {
        grid->data[i] = calloc(2 * grid->lh + COLS(data), sizeof(double));
        grid->data[i] += grid->lh;
    }
    grid->data += grid->lh;

    // Allocate G matrix with halo
    grid->G = calloc(2 * grid->lh + ROWS(data), sizeof(double*));
    for (i = 0; i < ROWS(data) + 2 * grid->lh; ++i) 
    {
        grid->G[i] = calloc(2 * grid->lh + COLS(data), sizeof(double));
        for (j = 0; j < 2 * grid->lh + COLS(data); ++j) grid->G[i][j] = 1.0;
        grid->G[i] += grid->lh;
    }
    grid->G += grid->lh;

    // Allocate u matrix with halo
    grid->u[0] = calloc(2 * grid->lh + ROWS(data), sizeof(double*));
    for (i = 0; i < ROWS(data) + 2 * grid->lh; ++i) 
    {
        grid->u[0][i] = calloc(2 * grid->lh + COLS(data) + 1, sizeof(double));
        grid->u[0][i] += grid->lh;
    }
    grid->u[0] += grid->lh;

    // Allocate v matrix with halo
    grid->u[1] = calloc(2 * grid->lh + ROWS(data) + 1, sizeof(double*));
    for (i = 0; i < ROWS(data) + 2 * grid->lh + 1; ++i) 
    {
        grid->u[1][i] = calloc(2 * grid->lh + COLS(data), sizeof(double));
        grid->u[1][i] += grid->lh;
    }
    grid->u[1] += grid->lh;

    // Copy Data from matrix to padded data object
    for (i = 0; i < ROWS(data); ++i)
    {
        memcpy(grid->data[i], data->data + i * COLS(data), COLS(data) * sizeof(double));
        if (G != NULL)
        {
            memcpy(grid->G[i], G->data + i * COLS(data), COLS(data) * sizeof(double));
        }
        else
        {
            for (j = 0; j < COLS(data); ++j) grid->G[i][j] = 1.0;
        }
        
        
    }


    // Calculate Courant Numbers u
    double** cu = grid->u[0];
    for (i = 0; i < ROWS(data); ++i)
    {
        double dx = x->data[1] - x->data[0];
        cu[i][0] = u->data[i*COLS(u)] * dt / dx;
        for (j = 0; j < COLS(data) - 1; ++j)
        {
            dx = x->data[j+1] - x->data[j];
            cu[i][j+1] = 0.5 * (u->data[i*COLS(u)+j+1] + u->data[i*COLS(u)+j]) * dt / dx;
        }
        dx = x->data[j] - x->data[j-1];
        cu[i][j+1] = u->data[i*COLS(u)+j] * dt / dx;
    }

    // Calculate Courant Numbers v
    double** cv = grid->u[1];
    for (j = 0; j < COLS(data); ++j)
    {
        double dy = y->data[1] - y->data[0];
        cv[0][j] = v->data[j] * dt / dy;
        for (i = 0; i < ROWS(data) - 1; ++i)
        {
            dy = y->data[i+1] - y->data[i];
            cv[i+1][j] = 0.5 * (v->data[(i+1)*COLS(v)+j] + v->data[i*COLS(v)+j]) * dt / dy;
        }
        dy = y->data[i] - y->data[i-1];
        cv[i+1][j] = v->data[i*COLS(v)+j] * dt / dy;
    }

    return grid;
}

void memcpy_arakawa_2d(mat* data, arakawa_2d* grid)
{
    assert(data->N == grid->N);
    assert(data->M == grid->M);
    unsigned long i;
    for (i = 0; i < ROWS(data); ++i)
        memcpy(data->data + i*COLS(data), grid->data[i], sizeof(double) * COLS(data));
}

mat* get_mat_arakawa(arakawa_2d* grid)
{
    mat* res = allocate_matNxM(grid->N, grid->M);
    memcpy_arakawa_2d(res, grid);
    return res;
}

mat* free_arakawa_2d(arakawa_2d* grid, const int get_field)
{
    mat* res = NULL;
    if (get_field) res = get_mat_arakawa(grid);
    long i;
    for (i = -(int)grid->lh; i < (int)grid->N + (int)grid->lh; ++i)
    {
        free(grid->data[i] - grid->lh);
        free(grid->u[0][i] - grid->lh);
        free(grid->u[1][i] - grid->lh);
        free(grid->G[i]    - grid->lh);
    }
    free(grid->u[1][i] - grid->lh);
    free(grid->data - grid->lh);
    free(grid->G    - grid->lh);
    free(grid->u[0] - grid->lh);
    free(grid->u[1] - grid->lh);
    free(grid);
    return res;
}

static double** empty_psi2d(arakawa_2d* grid)
{
    unsigned long i;
    double** res = calloc(grid->N + 2 * grid->lh, sizeof(double*));   
    for (i = 0; i < grid->N + 2 * grid->lh; ++i)
        res[i] = (double*)(calloc(grid->M + 2 * grid->lh, sizeof(double))) + grid->lh;

    return (res + grid->lh);
}

static double **copy_courant2d(arakawa_2d* grid, int axis)
{
    assert(axis == 0 || axis == 1);
    unsigned long i;
    if (axis == 1)
    {
        double** out = calloc(2 * grid->lh + grid->N, sizeof(double*));
        for (i = 0; i < grid->N + 2 * grid->lh; ++i)
        {
            out[i] = calloc(grid->M + 2 * grid->lh + 1, sizeof(double));
            memcpy(out[i], grid->u[0][i-grid->lh] - grid->lh, sizeof(double) * (2 * grid->lh + grid->M + 1));
            out[i] += grid->lh;
        }
        return (out + grid->lh);
    }

    double** out = calloc(2 * grid->lh + grid->N + 1, sizeof(double*));
    for (i = 0; i < grid->N + 2 * grid->lh + 1; ++i)
    {
        out[i] = calloc(grid->M + 2 * grid->lh, sizeof(double));
        memcpy(out[i], grid->u[1][i-grid->lh] - grid->lh, sizeof(double) * (2 * grid->lh + grid->M));
        out[i] += grid->lh;
    }
    return (out + grid->lh);
}

//#define FREE_PSI(X,O) {unsigned long __i; for (__i = -2; __i < (int)grid->N + 2 + O; ++__i) free((X)[__i]-2); free((X)-2);}
static void __free_2d_grid(double** data, int pad, unsigned NY)
{
    unsigned long i;
    data -= pad;
    for (i = 0; i < NY + 2 * pad; ++i) free(data[i] - pad);
    free(data);
}

arakawa_2d* mpadvec2d(arakawa_2d* grid, const long order, const unsigned options)
{
    assert(order >= 1);
    // Internal Variables
    long o, i, j;
    // FLAGS
    int FT = (options & MPDATA_NONOSCIL) ? 1 : 0;

    // Helpfull pointers
    double** psip = grid->data; 
    double** psin = empty_psi2d(grid);
    double** C[2] = {copy_courant2d(grid, 1), copy_courant2d(grid, 0)};
    double** G = grid->G;

    // Fill Halo
    assert(fill_halo_y(psip, C[1], G, grid->N, grid->M, options));
    for (i = 0; i < grid->N; ++i)
        fill_halo_x(psip[i], C[0][i], G[i], grid->M, options);

    // Temp
    double flow[2] = {0.0, 0.0};
    double** diffuv[2] = {copy_courant2d(grid, 1), copy_courant2d(grid, 0)};

    // non oscillatory data fields
    double** psi_max = NULL;
    double** psi_min = NULL;
    if (FT)
    {
        psi_max = empty_psi2d(grid);
        psi_min = empty_psi2d(grid);
        for (i = -1; i < grid->N + 1; ++i) for(j = -1; j < grid->M + 1; ++j)
        {
            psi_max[i][j] = va_max(5, psip[i-1][j], psip[i][j-1], psip[i][j], psip[i][j+1], psip[i+1][j]);
            psi_min[i][j] = va_min(5, psip[i-1][j], psip[i][j-1], psip[i][j], psip[i][j+1], psip[i+1][j]);
            //if (i < 1) printf("%ix%i\n", i, j);
        }

        // Fix the Lower and Upper bound
        for(j = -1; j < (int)grid->M + 1; ++j)
        {
            i = -2;
            psi_max[i][j] = va_max(4, psip[i][j-1], psip[i][j], psip[i][j+1], psip[i+1][j]);
            psi_min[i][j] = va_min(4, psip[i][j-1], psip[i][j], psip[i][j+1], psip[i+1][j]);
            i = grid->N+1;
            psi_max[i][j] = va_max(4, psip[i][j-1], psip[i][j], psip[i][j+1], psip[i-1][j]);
            psi_min[i][j] = va_min(4, psip[i][j-1], psip[i][j], psip[i][j+1], psip[i-1][j]);
        }

        // Fix Left and Right Bound
        for(i = -1; i < (int)grid->N + 1; ++i)
        {
            j = -2;
            psi_max[i][j] = va_max(4, psip[i-1][j], psip[i][j], psip[i][j+1], psip[i+1][j]);
            psi_min[i][j] = va_min(4, psip[i-1][j], psip[i][j], psip[i][j+1], psip[i+1][j]);
            j = grid->M+1;
            psi_max[i][j] = va_max(4, psip[i-1][j], psip[i][j-1], psip[i][j], psip[i+1][j]);
            psi_min[i][j] = va_min(4, psip[i-1][j], psip[i][j-1], psip[i][j], psip[i+1][j]);
        }

        // Edges
        psi_max[-2][-2] = va_max(3, psip[-2][-2], psip[-1][-2], psip[-2][-1]);              
        psi_min[-2][-2] = va_min(3, psip[-2][-2], psip[-1][-2], psip[-2][-1]); 

        j = grid->M+1;
        psi_max[-2][j] = va_max(3, psip[-2][j], psip[-1][j], psip[-2][j-1]);        
        psi_min[-2][j] = va_min(3, psip[-2][j], psip[-1][j], psip[-2][j-1]);     

        i = grid->N+1;
        psi_max[i][-2] = va_max(3, psip[i][-2], psip[i-1][-2], psip[i][-1]);         
        psi_min[i][-2] = va_min(3, psip[i][-2], psip[i-1][-2], psip[i][-1]);     

        psi_max[i][j] = va_max(3, psip[i][j], psip[i-1][j], psip[i][j-1]);         
        psi_min[i][j] = va_min(3, psip[i][j], psip[i-1][j], psip[i][j-1]); 

    }

    
    for (o = 0; o < order; ++o)
    {

        // Calculate cell Fluxes
        for (i = 0; i < (int)grid->N; ++i) for (j = 0; j < (int)grid->M; ++j)
        {
            flow[0]  = flux(psip[i  ][j], psip[i+1][j], C[1][i+1][j]);
            flow[0] -= flux(psip[i-1][j], psip[i  ][j], C[1][i  ][j]);
            flow[1]  = flux(psip[i][j  ], psip[i][j+1], C[0][i][j+1]);
            flow[1] -= flux(psip[i][j-1], psip[i][j  ], C[0][i][j  ]);
            psin[i][j] = psip[i][j] - (flow[0] + flow[1]) / G[i][j];
        }

        // Refill Halo 
        fill_halo_y(psin, C[1], G, grid->N, grid->M, options);
        for (i = 0; i < grid->N; ++i)
            fill_halo_x(psin[i], C[0][i], G[i], grid->M, options);

        // Determine Antidiffusive Velocites
        for (i = -1; i < (long)grid->N; ++i) for (j = -1; j < (int)grid->M; ++j)
        {
            double Aval, Bval;
            // V Velocities
            if (j >= 0)
            {
                Aval = mp_A(psin[i+1][j], psin[i][j]);
                Bval = mp_B2d(psin, C[1], C[0], i,j, 0);
                diffuv[1][i+1][j] = fabs(C[1][i+1][j]) * (1. - fabs(C[1][i+1][j])) * Aval - Bval;
            }

            // U Velocities
            if (i >= 0)
            {
                Aval = mp_A(psin[i][j+1], psin[i][j]);
                Bval = mp_B2d(psin, C[1], C[0], i,j, 1);
                diffuv[0][i][j+1] = fabs(C[0][i][j+1]) * (1. - fabs(C[0][i][j+1])) * Aval - Bval;
            }
        }

        if(!FT) goto skip_nonoscil;
        double** Cm[2] = {copy_courant2d(grid, 1), copy_courant2d(grid, 0)};
        for (i = -1; i < (long)grid->N; ++i) for (j = -1; j < (int)grid->M; ++j)
        {
            double betaf = 0.0;
            // V Velocities
            if (j < 0) goto skip_vcorr;
            if (diffuv[1][i+1][j] > 0.0)
                if (psin[i][j] > 0.0)
                    betaf = va_min(3, 1.0, beta_dn2d(psin, diffuv, psi_min, G, i, j), beta_up2d(psin, diffuv, psi_max, G, i+1, j));
                else
                    betaf = va_min(3, 1.0, beta_up2d(psin, diffuv, psi_max, G, i, j), beta_dn2d(psin, diffuv, psi_min, G, i+1, j));
            else
                if (psin[i+1][j] > 0.0)
                    betaf = va_min(3, 1.0, beta_up2d(psin, diffuv, psi_max, G, i, j), beta_dn2d(psin, diffuv, psi_min, G, i+1, j));
                else
                    betaf = va_min(3, 1.0, beta_dn2d(psin, diffuv, psi_min, G, i, j), beta_up2d(psin, diffuv, psi_max, G, i+1, j));
            Cm[1][i+1][j] = diffuv[1][i+1][j] * betaf;

            // U Velocities
            skip_vcorr:
            if (i < 0) continue;
            if (diffuv[0][i][j+1] > 0.0)
                if (psin[i][j] > 0.0)
                    betaf = va_min(3, 1.0, beta_dn2d(psin, diffuv, psi_min, G, i, j), beta_up2d(psin, diffuv, psi_max, G, i, j+1));
                else
                    betaf = va_min(3, 1.0, beta_up2d(psin, diffuv, psi_max, G, i, j), beta_dn2d(psin, diffuv, psi_min, G, i, j+1));
            else
                if (psin[i][j+1] > 0.0)
                    betaf = va_min(3, 1.0, beta_up2d(psin, diffuv, psi_max, G, i, j), beta_dn2d(psin, diffuv, psi_min, G, i, j+1));
                else
                    betaf = va_min(3, 1.0, beta_dn2d(psin, diffuv, psi_min, G, i, j), beta_up2d(psin, diffuv, psi_max, G, i, j+1));
            Cm[0][i][j+1] = diffuv[0][i][j+1] * betaf;
        }
        
        // Update Fields with Cm
        // Refill Halo 
        fill_halo_y(psin, Cm[1], G, grid->N, grid->M, options);
        for (i = 0; i < grid->N; ++i)
            fill_halo_x(psin[i], Cm[0][i], G[i], grid->M, options);
        // Update Fields
        for (i = -grid->lh; i < (int)(grid->N + grid->lh); ++i)
        {
            memcpy(psip[i], psin[i], sizeof(double) * grid->M);
            memcpy(C[0][i] - grid->lh, Cm[0][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh + 1));
            memcpy(C[1][i] - grid->lh, Cm[1][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh));
        }
        memcpy(C[1][i] - grid->lh, Cm[1][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh));

        __free_2d_grid(Cm[0], grid->lh, grid->N);
        __free_2d_grid(Cm[1], grid->lh, grid->N+1);
        continue;
        skip_nonoscil:

        // Update Fields
        fill_halo_y(psin, diffuv[1], G, grid->N, grid->M, options);
        for (i = 0; i < grid->N; ++i)
            fill_halo_x(psin[i], diffuv[0][i], G[i], grid->M, options);
        for (i = -grid->lh; i < (int)(grid->N + grid->lh); ++i)
        {
            memcpy(psip[i], psin[i], sizeof(double) * grid->M);
            memcpy(C[0][i] - grid->lh, diffuv[0][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh + 1));
            memcpy(C[1][i] - grid->lh, diffuv[1][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh));
        }
        memcpy(C[1][i] - grid->lh, diffuv[1][i] - grid->lh, sizeof(double) * (grid->M + 2 * grid->lh));

    }

    __free_2d_grid(psin,grid->lh,grid->N);
    __free_2d_grid(diffuv[0],grid->lh,grid->N);
    __free_2d_grid(diffuv[1],grid->lh,grid->N+1);
    __free_2d_grid(C[0],grid->lh,grid->N);
    __free_2d_grid(C[1],grid->lh,grid->N+1);
    if (FT)
    {
        __free_2d_grid(psi_max, grid->lh, grid->N);
        __free_2d_grid(psi_min, grid->lh, grid->N);
    }

    return grid;
}