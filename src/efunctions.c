#include "../include/environment.h"
#include "../include/extio.h"


const double ecoll[21][15] =
{
    {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
    0.001,0.001,0.001,0.001,0.001},{0.003,0.003,0.003,0.004,0.005,
    0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970},
    {0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430,
    0.580,0.790,0.930,0.960,1.000},{0.009,0.009,0.009,0.012,0.015,
    0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000},
    {0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770,
    0.840,0.950,0.970,1.000,1.000},{0.017,0.017,0.017,0.020,0.022,
    0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000},
    {0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870,
    0.900,0.950,1.000,1.000,1.000},{0.025,0.025,0.025,0.036,0.043,
    0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000},
    {0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900,
    0.940,1.000,1.000,1.000,1.000},{0.030,0.030,0.030,0.047,0.064,
    0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000},
    {0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910,
    0.950,1.000,1.000,1.000,1.000},{0.035,0.035,0.035,0.055,0.079,
    0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000},
    {0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910,
    0.950,1.000,1.000,1.000,1.000},{0.037,0.037,0.037,0.060,0.080,
    0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000},
    {0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920,
    0.950,1.000,1.000,1.000,1.000},{0.037,0.037,0.037,0.052,0.067,
    0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000},
    {0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950,
    1.000,1.000,1.000,1.000,1.000},{0.036,0.036,0.036,0.042,0.048,
    0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020},
    {0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030,
    1.040,1.040,1.040,1.040,1.040},{0.033,0.033,0.033,0.033,0.033,
    0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300},
    {0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000,
    4.000,4.000,4.000,4.000,4.000}
};

vec* convert_z_to_p(vec* z, vec* T, double p0)
{
    assert(z->N == T->N);
    assert(p0 > 0.0);
    vec* p = fill_vec(z->N, p0);
    unsigned long i;
    // const double scale_H = -GVAL / RAIR;
    for (i = 1; i < T->N; ++i)
    {
        vec* Tinv = slice_vec(T, 0, i);
        vec* zslc = slice_vec(z, 0, i);

        if (LEN(Tinv) > 1)
        {
            //pow_vec(Tinv, -1.0);
            //p->data[i] = p0 * exp(scale_H * trapz(Tinv, zslc));
            double dz = z->data[i] - z->data[i-1];
            p->data[i] = p->data[i-1] - p->data[i-1] * GVAL / RAIR / T->data[i-1] * dz;
        }
        else
        {
            double dz = z->data[i] - z->data[0];
            p->data[i] = p0 - p0 * GVAL / RAIR / T->data[i-1] * dz;
        }
        

        free_vec(Tinv);
        free_vec(zslc);
    }


    return p;
}

vec* convert_p_to_z(vec* p, vec* T, double z0)
{
    assert(p->N == T->N);
    const double scale_H = -RAIR / GVAL;
    double dp = 0.0;
    vec* z = fill_vec(p->N, z0);
    unsigned long i;

    for (i = 1; i < p->N; ++i)
    {
        dp = p->data[i] - p->data[i-1];
        z->data[i] = z->data[i-1] + scale_H / p->data[i] * T->data[i] * dp;
    }

    return z;
}

/*Enviroment Calculations*/
mat* ideal_gass_pressure(mat* T0, vec* z0, const double p0)
{
    mat* po0 = fill_mat(ROWS(T0), COLS(T0), p0);
    const double scaleH = -GVAL / RAIR;
    double integral;
    unsigned long i, j;
    unsigned long ZLEN = LEN(z0);

    for (i = 0; i < COLS(po0); ++i)
    {
        vec* _tmp = slice_mat_col(T0, i);
        pow_vec(_tmp, -1.0);

        for (_tmp->N = 1, j = 1; j < ROWS(T0); ++j, ++_tmp->N)
        {
            z0->N = j;
            integral = p0 * exp(scaleH * trapz(_tmp, z0));
            set_mat_ind(po0, j, i, integral);
        }

        free_vec(_tmp);
    }
    z0->N = ZLEN;

    return po0;
}

mat* ideal_gass_density(mat* T0, vec* z0, const double p0)
{
    mat* rho0 = fill_mat(ROWS(T0), COLS(T0), p0 / T0->data[0] / RAIR);
    const double scaleH = -GVAL / RAIR;
    double integral;
    unsigned long i, j;
    unsigned long ZLEN = LEN(z0);
    
    for (i = 0; i < COLS(rho0); ++i)
    {
        vec* _tmp = slice_mat_col(T0, i);
        pow_vec(_tmp, -1.0);

        for (_tmp->N = 1, j = 1; j < ROWS(T0); ++j, ++_tmp->N)
        {
            z0->N = j;
            integral = p0 * exp(scaleH * trapz(_tmp, z0));
            set_mat_ind(rho0, j, i, integral / RAIR / T0->data[j * COLS(T0) + i]);
        }

        free_vec(_tmp);
    }
    z0->N = ZLEN;

    return rho0;
}

mat* potential_T_converter(mat* rho, mat* T)
{
    mat* cfnc = copy_mat(rho);
    unsigned long i;
    for (i = 0; i < LEN_MAT(rho); ++i)
        cfnc->data[i] *= RAIR * T->data[i];
    
    double p0 = max_mat(cfnc);
    for (i = 0; i < LEN_MAT(rho); ++i)
        cfnc->data[i] = pow(p0 / cfnc->data[i], RAIR / CPA);
    return cfnc;
}

/*Kinematic Framework*/
void strmf_shallow_cumulus(vec* x, vec* y, double t, mat* strm)
{
    assert(y->N == strm->N); assert(x->N == strm->M);
    assert(t >= 0.0);
    const double hx = 1.8e3;
    const double xc = 4.5e3;
    const double x0 = 3.6e3;
    const double Lz = 2.7e3;

    unsigned long i,j;
    t = MIN(t, 2400.0);
    memset(strm->data, 0, sizeof(double) * LEN_MAT(strm));

    // Determine Coefficients 
    double A1 = 5.73e2, A2 = 0.0;
    double alpha, beta;
    double z0, hz;
    double T1a, T1b, T2;
    double res;

    if (t > 300.0 && t <= 900.0)
    {
        A2 = 6e2 * (1 + cos(PI * ((t - 300.0) / 600.0 - 1.)));
    }
    else if (t > 900.0 && t <= 1500.0)
    {
        A1 = 5.73e2 + 2.02e3 * (1 + cos(PI * ((t - 900.0) / 600.0 + 1)));
        A2 = 6e2 * (1 + cos(PI * (t - 300.0) / 600.0 - 1.));
    }
    else if (t > 1500.0)
    {
        A1 = 1.15e3 + 1.72e3 * (1 + cos(PI * ((t - 1500.0) / 900.0)));
        A2 = 5.00e2 * (1 + cos(PI * ((t - 1500.0) / 900.0 - 1.)));
    }

    for (i = 0; i < LEN(y); ++i)
    {
        if (y->data[i] >= 2.7e3) continue;
        z0 = (y->data[i] <= 1.7e3) ? 0.0 : 7.2e2;
        hz = (y->data[i] <= 1.7e3) ? 3.4e3 : 2e3;

        for (j = 0; j < LEN(x); ++j)
        {
            alpha = (fabs(x->data[j] - xc) <= 0.9e3) ? 1.0 : 0.0;
            beta = (x->data[j] <= 5.4e3) ? 1.0 : -1.0;

            T1a = cos(alpha * PI * (x->data[j] - x0) / hx);
            T1b = sin(beta * PI * (y->data[i] - z0) / hz);
            T2  = (y->data[i] / Lz) * (y->data[i] / Lz);

            res = -A1 * T1a * T1b + A2 / 2.0 * T2;
            set_mat_ind(strm, i, j, res);
        }
    }
    //fprint_vec(stdout, y);
}

void strmf_stratocumulus_sym(vec* x, vec* y, double t, mat* strm)
{
    assert(LEN(x) == COLS(strm));
    assert(LEN(y) == ROWS(strm));
    unsigned long i, j;
    double T1, T2;
    const double A = MODEL_SETTINGS.strm_density;
    const double Ztop = MODEL_SETTINGS.Ztop;
    const double Xwth = MODEL_SETTINGS.Xwidth;
    for (i = 0; i < ROWS(strm); ++i) for (j = 0; j < COLS(strm); ++j)
    {
        T2 = cos(2.0 * PI * x->data[j] / Xwth);
        T1 = sin(PI * y->data[i] / Ztop);
        set_mat_ind(strm, i, j, -A * T1 * T2);
        
    }
}

void strmf_stratocumulus_asym(vec* x, vec* y, double t, mat* strm)
{
    assert(LEN(x) == COLS(strm));
    assert(LEN(y) == ROWS(strm));
    unsigned long i, j;
    double T1, T2;
    const double A = MODEL_SETTINGS.strm_density;
    const double Ztop = MODEL_SETTINGS.Ztop;
    const double Zclb = MODEL_SETTINGS.Zclb;
    const double Xwth = MODEL_SETTINGS.Xwidth;

    for (i = 0; i < ROWS(strm); ++i) for (j = 0; j < COLS(strm); ++j)
    {
        T2 = cos(2.0 * PI * x->data[j] / Xwth);
        if (y->data[i] <= Zclb)
            T1 = sin(PI / 2.0 * y->data[i] / Zclb);
        else
            T1 = sin(PI / 2.0 * ((y->data[i] - Zclb) / (Ztop - Zclb) + 1.0));
        set_mat_ind(strm, i, j, -A * T1 * T2);
        
    }
}

/*Droplet Physics*/
double droplet_vT(double r)
{
    double dmu = r * 2e6;
    double alpha = 4.5795e5, beta = 2. / 3.;
    double vT;

    if (dmu >= 134.43 && dmu < 1511.64)
    {
        alpha = 4.962e3;
        beta = 1. / 3.;
    }
    else if (dmu >= 1511.64 && dmu < 3477.84)
    {
        alpha = 1.732e3;
        beta = 1. / 6.;
    }
    else if (dmu >= 3477.84)
    {
        alpha = 9.17e2;
        beta = 0.0;
    }
    
    vT = alpha * pow((droplet_mass(r)) * 1e3, beta); // Is in cm/s
    return vT * 1e-2;

}

double ventilation_coef(double r)
{
    double vT = droplet_vT(r);
    const double Nsh = 0.71;
    double Nre = 2.0 * r * vT / KIN_AIR;
    double fac = pow(Nsh, 1. / 3.) * sqrt(Nre);
    if (fac < 1.4)
        return 1.00 + 0.108 * fac * fac;
    return 0.78 + 0.308 * fac;
}

vec* droplet_growth_rate(vec* r)
{
    unsigned long i;
    vec* dr_dt = zero_vec(LEN(r));
    const double Av = 1e-10;

    for (i = 0; i < LEN(r); ++i)
        dr_dt->data[i] = Av * ventilation_coef(r->data[i]) / r->data[i];

    return dr_dt;
}

void courant_coal(vec* r, mat* c, int* ima)
{
    assert(LEN_MAT(c) == LEN(r) * LEN(r));
    unsigned long i, j, k;
    unsigned long N = LEN(r);
    int kk;
    double x0, xi, xj, xk, xkm;

    memset(c->data, 0, sizeof(double) * N * N);
    memset(ima, 0, sizeof(int) * N * N);

    //mat* edg= edges(r);
    vec* dlnr = zero_vec(N);
    for (i = 0; i < N-1; ++i)
    {
        dlnr->data[i] = log(r->data[i+1] / r->data[i]);
    }
    dlnr->data[i] = dlnr->data[i-1];

    for (i = 0; i < N; ++i)
    {
        for (j = i; j < N; ++j)
        {
            xi = droplet_mass(r->data[i]);
            xj = droplet_mass(r->data[j]);
            x0 = xi + xj;
            for (k = j; k < N; ++k)
            {
                kk = -1;
                xk = droplet_mass(r->data[k]);
                xkm = (k > 0) ? droplet_mass(r->data[k-1]) : droplet_mass(r->data[N-1]);
                if (xk >= x0 && xkm < x0)
                {
                    if (c->data[i * N + j] < 1.0 - 1e-8)
                    {
                        kk = (int)k - 1;
                        c->data[i * N + j] = log(x0 / xkm) / (3.0 * dlnr->data[i]);
                    }
                    else
                    {
                        c->data[i * N + j] = 0.0;
                        kk = (int)k;
                    }
                    ima[i*N + j] = MIN((int)N - 2, kk);
                    break;                   
                }
            }
            c->data[j*N+i] = c->data[i*N+j];
            ima[j*N + i] = ima[i*N + j];
        }
    }
    free_vec(dlnr);
}

mat* collision_kernel(vec* r, const int type, const double dt)
{
    unsigned long i, j;
    unsigned long N = LEN(r);
    double xi, xj, E;
    mat* ck = zero_mat(N,N);

    if (type == KERNEL_GOLV)
    {
        for (i = 0; i < N; ++i) for (j = 0; j < i; ++j)
        {
            xi = droplet_mass(r->data[i]);
            xj = droplet_mass(r->data[j]);
            ck->data[i * N + j] = 1.5 * (xi + xj);
            ck->data[j * N + i] = ck->data[i * N + j];
        }
    }
    else if (type == KERNEL_LONG)
    {
        for (i = 0; i < N; ++i) for (j = 0; j < i; ++j)
        {
            xi = droplet_mass(r->data[i]);
            xj = droplet_mass(r->data[j]);
            E = 1.0;
            if (r->data[i] < 50e-6)
            {
                double rx = r->data[i] * 1e2;
                double ry = r->data[j] * 1e2;
                E = MAX(4.5e4 * rx * rx * (1. - 3e-4 / ry), 1e-3);
            }

            ck->data[i * N + j] = PI * (r->data[i] + r->data[j]) * (r->data[i] + r->data[j]) *\
                                 E * fabs(droplet_vT(r->data[i]) - droplet_vT(r->data[j]));
            ck->data[j * N + i] = ck->data[i * N + j];
        }
    }
    else if (type == KERNEL_HALL)
    {
        const double r0[] = {6e-6, 8e-6, 10e-6, 15e-6, 20e-6, 25e-6, 30e-6,
                            40e-6, 50e-6, 60e-6, 70e-6, 100e-6, 150e-6, 200e-6, 300e-6};
        const double Ra[] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.50,
                            0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
        unsigned long k, ir, iq;

        for (i = 0; i < N; ++i) for(j = 0; j < i; ++j)
        {
            iq = 0; ir = 0;
            for (k = 1; k < 15; ++k)
            {
                if (r->data[j] <= r0[k] && r->data[j] >= r0[k-1])
                    ir = k;
                else if (r->data[j] > r0[14])
                    ir = 15;
                else if (r->data[j] < r0[0])
                    ir = 1;
            }

            double cra = r->data[j] / r->data[i];
            for (k = 1; k < 21; ++k)
                if (cra < Ra[k] && cra > Ra[k-1]) iq = k;

            if (ir < 15)
            {
                if (ir >= 1)
                {
                    double p = (r->data[j] - r0[ir-1]) / (r0[ir] - r0[ir-1]);
                    double q = (cra - Ra[iq-1]) / (Ra[iq] - Ra[iq-1]);
                    ck->data[j * N + i] = (1.-p)*(1.-q)*ecoll[iq-1][ir-1] +\
                        p*(1.-q)*ecoll[iq-1][ir] +\
                        q*(1.-p)*ecoll[iq][ir-1] +\
                        p*q*ecoll[iq][ir];

                }
                else
                {
                    double q = (cra - Ra[iq-1]) / (Ra[iq] - Ra[iq-1]);
                    ck->data[j * N + i] = MIN((1.-q)*ecoll[iq-1][0] + q*ecoll[iq][0], 1.0);
                }
                
            }
            else
            {
                double q = (cra - Ra[iq-1]) / (Ra[iq] - Ra[iq-1]);
                ck->data[j * N + i] = MIN((1.-q)*ecoll[iq-1][14] + q*ecoll[iq][14], 1.0);
            }
            ck->data[i * N + j] = ck->data[j * N + i];
            
        }

        for (i = 0; i < N; ++i) for(j = 0; j < i; ++j)
        {
            ck->data[j*N+i] = PI * (r->data[i] + r->data[j]) * (r->data[i] + r->data[j]) * ck->data[j*N+i] *\
                             fabs(droplet_vT(r->data[j]) - droplet_vT(r->data[i]));
            ck->data[i*N+j] = ck->data[j*N+i];
        }
        
    }
    
    // Smooth Collision Kernel
    int jm, im, jp, ip;
    double dlnr;
    mat* cko = allocate_matNxN(N);

    for (i = 0; i < N; ++i)
    {
        im = MAX((int)i - 1, 0);
        ip = MIN((int)i + 1, (int)N - 1);
        dlnr = (i + 1 == N) ? log(r->data[i] / r->data[i-1]) : log(r->data[i+1] / r->data[i]);

        for (j = 0; j < N; ++j)
        {
            jm = MAX((int)j - 1, 0);
            jp = MIN((int)j + 1, (int)N - 1);
        

            cko->data[i * N + j] = 0.125 * (ck->data[im * N + j] + ck->data[i * N + jm] + ck->data[ip * N + j] + ck->data[i * N + jp]);
            cko->data[i * N + j] += 0.5 * ck->data[i * N + j];
            if (i == j) cko->data[j * N + i] *= 0.5;

            cko->data[i * N + j] *= dt * dlnr;
        }
    }

    free_mat(ck);
    return cko;
}

vec* collision(vec* bin, vec* r, int* ima, mat* c, mat* ck)
{
    long N = (long)LEN(r);
    long j, i, i0, il, k, kp;
    const double gmin = 1e-60;
    double xi, xj, x0;
    double gsi, gsj, gsk, gk;

    i0 = 0; il = N - 1;

    // Needed For Conversion
    mat* bedge = edges(r);
    vec* dlnr = allocate_vec(N);
    vec* dr = dr_bins(r);


    for (i = 0; i < N; ++i)
        dlnr->data[i] = log(bedge->data[2 * i + 1] / bedge->data[2 * i]);
    free_mat(bedge);

    // Convert to mass density
    for (i = 0; i < N; ++i) 
        bin->data[i] *= droplet_mass(r->data[i]) * dr->data[i] / dlnr->data[i];
    
    // Collision Routine
    // * Find Upper and Lower bound
    for (i = 0; i < N - 1; ++i)
    {
        i0 = i;
        if (bin->data[i] > gmin) break;
    }

    for (i = N-2; i > 0; --i)
    {
        il = i;
        if (bin->data[i] > gmin) break;
    }

    // Do Droplet Interaction
    for (i = i0; i < il; ++i)
        {
        for (j = i; j < il; ++j)
        {
            k = ima[i * N + j];
            kp = k + 1;
            xi = droplet_mass(r->data[i]);
            xj = droplet_mass(r->data[j]);

            x0 = ck->data[i * N + j] * bin->data[i] * bin->data[j];
            x0 = MIN(x0, bin->data[i] * xj);
            if (j != k) x0 = MIN(x0, bin->data[j] * xi);

            gsi = x0 / xj;
            gsj = x0 / xi;
            gsk = gsi + gsj;
            bin->data[i] = bin->data[i] - gsi;
            bin->data[j] = bin->data[j] - gsj;
            gk = bin->data[k] + gsk;
            
            if (gk > gmin)
            {
                double xl = log(bin->data[kp] / gk + 1e-60);
                double flx = gsk / xl* (exp(0.5*xl) - exp(xl*(0.5 - c->data[i * N + j])));
                flx = MIN(flx, gk);
                bin->data[k] = gk - flx;
                bin->data[kp] = bin->data[kp] + flx;
            }
        }
    }

    // Convert back to number  density
    for (i = 0; i < N; ++i) 
        bin->data[i] *= dlnr->data[i] / droplet_mass(r->data[i]) / dr->data[i];
    free_vec(dlnr); free_vec(dr);
    return bin;
}

/*Bin Definition and Conversion Formulas*/
vec* linexp_grid(const unsigned long N, const double linT, const double expT)
{
    vec* v = allocate_vec(N);
    unsigned long i; 
    for (i = 0; i < N; ++i)
        v->data[i] = (linT * (double)i + pow(10.0, (double)i * expT)) * 1e-6;
    return v;
}

vec* linmass_grid(const unsigned long N, const double linT, const double massS)
{
    vec* v = allocate_vec(N);
    double mi = droplet_mass(1e-6);
    const double mf = pow(2., 1. / massS);

    unsigned long i;
    for (i =  0; i < N; ++i)
    {
        v->data[i] = linT * (double)i * 1e-6;
        mi *= mf;
        v->data[i] += pow(3. / 4. / PI / RHOW * mi, 1. / 3.);
    }

    return v;
}

mat* edges(const vec* r)
{
    mat* edg = zero_mat(LEN(r), 2);
    unsigned long i;
    double ledge, redge;

    set_mat_ind(edg, 0, 1, 0.5 * (r->data[1] + r->data[0]));
    ledge = MAX(r->data[0] - (edg->data[1] - r->data[0]), 1e-10);
    set_mat_ind(edg, 0, 0, ledge);

    //printf("le %.2e re %.2e\n",  ledge, 0.5 * (r->data[1] + r->data[0]));

    for (i = 1; i < LEN(r) - 1; ++i)
    {
        redge = 0.5 * (r->data[i+1] + r->data[i]);
        ledge = 0.5 * (r->data[i] + r->data[i-1]);

        set_mat_ind(edg, i, 0, ledge);
        set_mat_ind(edg, i, 1, redge);
    }

    set_mat_ind(edg, i, 0, 0.5 * (r->data[i] + r->data[i-1]));
    redge = 2.0 * r->data[LEN(r) - 1] - edg->data[2 * (LEN(r) - 1)];
    set_mat_ind(edg, i, 1, redge);

    //printf("le %.2e re %.2e\n",  0.5 * (r->data[i] + r->data[i-1]), redge);

    return edg;
}

vec* dr_bins(const vec* r)
{
    vec* dr = zero_vec(LEN(r));
    mat* edg = edges(r);
    unsigned long i;
    for (i = 0; i < LEN(r); ++i)
        dr->data[i] = edg->data[2 * i + 1] - edg->data[2 * i];

    free_mat(edg);
    return dr;
}

/*Subscale Turbulence*/
vec* horizontal_diffusion(vec* v, vec* rho, vec* x, const double dt, int cyclic_field)
{
    assert(LEN(v) == LEN(rho));
    assert(LEN(v) == LEN(x));

    vec* turb = zero_vec(LEN(v));
    double Kscaler = 0.2 * sqrt(TRBE);
    unsigned long i;
    unsigned long N = LEN(v);
    int constant = 1;
    double dx = x->data[1] - x->data[0];
    double rhoc = rho->data[0];
    double K;
    
    for (i = 1; i < N; ++i)
        if (dx != (x->data[i] - x->data[i-1]) || rhoc != rho->data[i])
        {
            constant = 0;
            break;
        }

    if (constant)
    {
        K = Kscaler * dx;
        if (cyclic_field)
        {
            turb->data[0] = K * (v->data[v->N-1] + v->data[1] - 2.0 * v->data[0]) / dx / dx * dt;
            for (i = 1; i < N - 1; ++i)
                turb->data[i] = K * (v->data[i+1] + v->data[i-1] - 2.0 * v->data[i]) / dx / dx * dt;
            turb->data[i] = K * (v->data[i-1] + v->data[0] - 2.0 * v->data[i]) / dx / dx * dt;
        }
        else
        {
            turb->data[0] = K * (v->data[2] + v->data[0] - 2.0 * v->data[1]) / dx / dx * dt;
            for (i = 1; i < N - 1; ++i)
                turb->data[i] = K * (v->data[i+1] + v->data[i-1] - 2.0 * v->data[i]) / dx / dx * dt;
            turb->data[i] = K * (v->data[i] + v->data[i-2] - 2.0 * v->data[i-1]) / dx / dx * dt;
        }     
    }
    else
    {

        // We need to checl this further
        vec* vgrad = gradient_vec(v, x);
        mlt_vec(vgrad, 1, rho);
        for (i = 1; i < N; ++i)
            turb->data[i] = Kscaler * (x->data[i] - x->data[i-1]);
        turb->data[0] = Kscaler * (x->data[1] - x->data[0]);
        mlt_vec(turb, 1, vgrad);

        free_vec(vgrad);
        vgrad = gradient_vec(turb, x);
        free_vec(turb);
        turb = vgrad;
    }
    

    add_vec(v, 1, turb);
    free_vec(turb);
    return v;
}

/*Surface Parametrisation Fluxes*/
mat* SHF_T_flux(mat* T, const mat* rho, const double Hs, const double dt)
{
    double SHF = -MODEL_SETTINGS.LHF;
    double res;
    unsigned long i, j;

    for (j = 0; j < ROWS(T); ++j)
        for (i = 0; i < COLS(T); ++i)
        {
            res = -SHF / CPA / rho->data[j * COLS(rho) + i] / Hs * dt;
            T->data[j * COLS(T) + i] += res;
            T->data[j * COLS(T) + i] = MAX(0.0, T->data[j * COLS(T) + i]);
        }

    return T;
}

mat* LHF_q_flux(mat* q, const mat* rho, const double Hs, const double dt)
{
    double res;
    unsigned long i, j;

    for (j = 0; j < ROWS(q); ++j)
        for (i = 0; i < COLS(q); ++i)
        {
            res = -MODEL_SETTINGS.LHF / LV / rho->data[j * COLS(rho) + i] / Hs * dt;
            q->data[j * COLS(q) + i] += res;
            q->data[j * COLS(q) + i] = MAX(0.0, q->data[j * COLS(q) + i]);
        }

    return q;
}
