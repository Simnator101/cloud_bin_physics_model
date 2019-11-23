#include "../include/tests.h"
#include "../include/extio.h"

int test_vector_funcs()
{
    vec* tx = range(0, 12, 2);
    vec* ty = range(0, 6, 1);
    int rval = 0;

    printf("### Vector Functionality Tests ###\n# Vector data in x and y\n");
    fprint_vec(stdout, tx); 
    fprint_vec(stdout, ty);

    // Trapezium Test
    fprintf(stdout, "\nTrapezium Rule Result %f\n", trapz(ty, tx));
    if (trapz(ty, tx) != 25.0)
    {
        printf("Trapz invalid value of %.1f, should be 25.0\n", trapz(ty, tx));
        rval = 1;
    }

    // Gradient Test
    printf("\nGradient vector test dy/dx\n");
    vec* dtydtx = gradient_vec(ty, tx);
    fprint_vec(stdout, dtydtx);
    if (sum_vec(dtydtx) != 3.0)
    {
        printf("Sum of gradient values should be 3.0, but is %.1f\n", sum_vec(dtydtx));
        rval = 2;
    }

    free_vec(dtydtx);

    free_vec(tx);
    free_vec(ty);

    return rval;
}

int test_matrix_funcs()
{
    return 0;
}

int test_sounding()
{
    // Model data fields
    sounding snd = {NULL, NULL, NULL, NULL};

    // Initialise Sounding
    unsigned long long snd_id = (SOUNDING_Q << 16) | (SOUNDING_T << 8) | SOUNDING_P;
    read_sounding_txt(&snd, "./data/test_profile.txt", snd_id);
    printf("#### Sounding Info P(Pa), Z(m), T(K), Tde(K) ####\n");
    fprint_vec(stdout, snd.pe);
    fprint_vec(stdout, snd.ze);
    fprint_vec(stdout, snd.Te);
    fprint_vec(stdout, snd.Tde);
    

    double TeA[] = {298.350000, 297.950000, 296.750000, 295.650000, 294.950000, 293.650000,
                    293.050000, 291.350000, 289.950000, 287.950000, 286.450000,
                    285.050000, 284.150000, 284.450000, 284.050000, 284.350000,
                    283.350000, 284.150000, 284.350000, 283.150000, 281.950000,
                    279.950000, 275.950000, 269.950000};   
    if (memcmp(snd.Te->data, TeA, sizeof(double) * LEN(snd.Te)) != 0)
    {
        printf("WARNING: Test Sounding was incorrectly loaded.\n");
        free_sounding(&snd);
        return 1;
    }

    free_sounding(&snd);
    return 0;
}

void test_advection()
{
    FILE *pf = fopen("./tests/advection_test.txt", "w");
    unsigned long i;
    vec* v = fill_vec(100, 0.5);
    vec* x = range(0, (double)v->N, 1); 
    vec* c = zero_vec(100);
    for (i = 0; i < 100; ++i)
        c->data[i] = ((x->data[i] > 40.0 && x->data[i] <= 60.0) ? 1.0 : 0.0);

    // Print Out Original Concentration
    for (i = 0; i < LEN(c); ++i) fprintf(pf, "%.6f ", c->data[i]);
    fprintf(pf, "\r\n");

    // Run Simulation
    arakawa_1d* grid = make_arakawa_1d(c, x, v, NULL, 1.0);
    for (i = 0; i < (int)(100.0 / 0.5); ++i) mpadvec1d(grid, 2, MPDATA_CYCLIC_X | MPDATA_NONOSCIL);
    memcpy_arakawa_1d(c, grid);
    free_arakawa_1d(grid, 0);

    // Print Out Final Concentration
    for (i = 0; i < LEN(c); ++i) fprintf(pf, "%.6f ", c->data[i]);
    fprintf(pf, "\r\n");

    free_vec(c);
    free_vec(v);
    free_vec(x);
    fclose(pf);
}

void test_advection_2d()
{
    vec* x = linspace(0.0, 100.0, 101);
    vec* y = linspace(0.0, 100.0, 101);
    mat* u = allocate_matNxN(101);
    mat* v = allocate_matNxN(101);
    mat* c = allocate_matNxN(101);
    unsigned long i,j;
    FILE* pf = fopen("./tests/test_2dadv.txt", "w");

    for (i = 0; i < 101; ++i)
        for (j = 0; j < 101; ++j)
        {
            set_mat_ind(u, i, j,  0.1 * (y->data[i] - 50.0));
            set_mat_ind(v, i, j, -0.1 * (x->data[j] - 50.0));
            //set_mat_ind(u, i, j,  0.1);
            //set_mat_ind(v, i, j,  0.0 * (y->data[i] - 50.0));
            double T1 = pow((x->data[j] - 75), 2.0) / 200.0;
            double T2 = pow((y->data[i] - 50), 2.0) / 200.0;
            set_mat_ind(c, i, j, 3.0 * exp(-T1) * exp(-T2));
        }

    // Print Initial Concentration
    for (i = 0; i < 101 * 101; ++i)
        fprintf(pf, "%.6f ", c->data[i]);
    fprintf(pf, "\r\n");

    arakawa_2d* grid = make_arakawa_2d(c, x, y, u, v, NULL, 0.05);
    unsigned long NT = (unsigned long)(20.0 * PI / 0.05 * 12.0);
    for (i = 0; i < NT; ++i) mpadvec2d(grid, 4, MPDATA_CYCLIC_X | MPDATA_CYCLIC_Y | MPDATA_NONOSCIL);
    memcpy_arakawa_2d(c, grid);
    free_arakawa_2d(grid, 0);

    // Print End Concentration
    for (i = 0; i < 101 * 101; ++i)
        fprintf(pf, "%.6f ", c->data[i]);
    fprintf(pf, "\r\n");
    fclose(pf);

    free_vec(x);
    free_vec(y);
    free_mat(u);
    free_mat(v);
    free_mat(c);
}

void test_droplet_advc()
{
    unsigned long i, j;
    const double L = 1e-3;
    const double mm = droplet_mass(10e-6);
    const double pf = L / mm / mm;
    double T = 295.0;
    double p = 1e5;
    double dt = 0.5;
    unsigned long N = (unsigned long)(20.0 / dt);

    vec* r = linexp_grid(69, 0.25, 0.055);
    vec* dr = dr_bins(r);
    mat* edg = edges(r);
    mat* ba = zero_mat(N, LEN(r));
    vec* drdt = droplet_growth_rate(r);
    vec* q = fill_vec(N, 0.02);

    for (i = 0; i < LEN(r); ++i)
    {
        double mc = droplet_mass(r->data[i]);
        double dln = log(edg->data[2*i + 1] / edg->data[2*i]);
        ba->data[i] = 3.0 * mc * pf * exp(-mc / mm) * dln / dr->data[i];
    }
    free_mat(edg);

    FILE *fl = fopen("./tests/test_dadvec.txt", "w");
    fprintf(fl, "%.8f ", q->data[0]);

    for (i = 1; i < N; ++i)
    {
        if (i > N / 2) T = 298.15;
        double S = q->data[i-1] / sat_mixr_vapour(p, T) - 1.0;
        vec* b = slice_mat_row(ba, i-1);
        vec* C = courant_bin(drdt, dr, S, dt);

        vec* bn = advec_bin(b, C);
        for (j = 0; j < LEN(b); ++j) C->data[j] = (bn->data[j] - b->data[j]) / dt * droplet_mass(r->data[j]);

        q->data[i] = q->data[i-1] - trapz(C, r) * dt;
        set_mat_row(ba, bn, i);
        fprintf(fl, "%.8f ", q->data[i]);

        free_vec(C);
        free_vec(bn);
        free_vec(b);
    }

    
    fclose(fl);


    free_vec(dr);
    free_vec(r);
    free_vec(drdt);
    free_mat(ba);
    free_vec(q);
}

void test_droplet_coal()
{
    unsigned long i, j;
    const double L = 1e-3;
    const double mm = droplet_mass(10e-6);
    const double pf = L / mm / mm;
    vec* r = linexp_grid(120, 0.125, 0.032);
    vec* dr = dr_bins(r);
    mat* edg = edges(r);
    vec* b = zero_vec(LEN(r));
    fprint_vec(stdout, r);

    FILE* fp = fopen("./tests/test_golv.txt", "w");
    for (i = 0; i < LEN(r); ++i)
    {
        double mc = droplet_mass(r->data[i]);
        double dln = log(edg->data[2*i + 1] / edg->data[2*i]);
        b->data[i] = 3.0 * mc * pf * exp(-mc / mm) / dr->data[i] * dln;
    }

    int* ima = malloc(sizeof(int) * LEN(r) * LEN(r));
    mat* C = allocate_matNxN(LEN(r));
    courant_coal(r, C, ima);
    mat* cck = collision_kernel(r, KERNEL_LONG, 1.0);

    // Original dist
    for (i = 0; i < LEN(r); ++i)
    {
        double dln = log(edg->data[2*i + 1] / edg->data[2*i]);
        double mc = droplet_mass(r->data[i]);
        fprintf(fp, "%f ", b->data[i] / dln * dr->data[i] * mc);
    }
    fprintf(fp, "\r\n");

    for (i = 0; i < 1800; ++i)
        b = collision(b, r, ima, C, cck);

    // +30 Min dist
    for (i = 0; i < LEN(r); ++i)
    {
        double dln = log(edg->data[2*i + 1] / edg->data[2*i]);
        double mc = droplet_mass(r->data[i]);
        fprintf(fp, "%f ", b->data[i] / dln * dr->data[i] * mc);
    }
    fprintf(fp, "\r\n");

    // +60 min dist
    for (i = 0; i < 1800; ++i)
        b = collision(b, r, ima, C, cck);

    for (i = 0; i < LEN(r); ++i)
    {
        double dln = log(edg->data[2*i + 1] / edg->data[2*i]);
        double mc = droplet_mass(r->data[i]);
        fprintf(fp, "%f ", b->data[i] / dln * dr->data[i] * mc);
    }
    fclose(fp);

    // Print Out IMA CCK info
    fp = fopen("./tests/test_ima.txt", "w");
    for (i = 0; i < LEN(r); ++i) 
    {
        for (j = 0; j < LEN(r); ++j)
            fprintf(fp, "%i ", ima[i * LEN(r) + j]);
        fprintf(fp, "\r\n");
    }
    fclose(fp);
    
    fp = fopen("./tests/test_colC.txt", "w");
    for (i = 0; i < LEN(r); ++i) 
    {
        for (j = 0; j < LEN(r); ++j)
            fprintf(fp, "%.6e ", C->data[i * LEN(r) + j]);
        fprintf(fp, "\r\n");
    }
    fclose(fp);

    fp = fopen("./tests/test_ckern.txt", "w");
    for (i = 0; i < LEN(r); ++i) 
    {
        for (j = 0; j < LEN(r); ++j)
            fprintf(fp, "%.6e ", cck->data[i * LEN(r) + j]);
        fprintf(fp, "\r\n");
    }
    fclose(fp);

    free(ima);
    free_mat(cck);
    free_mat(C);
    free_vec(r);
    free_vec(dr);
    free_mat(edg);
    free_vec(b);
}

void test_eddy_mixing()
{
    unsigned long i, j;
    mat* C = allocate_matNxN(100);
    vec* z = range(0.0, 1000.0, 10.0);
    vec* x = range(0.0, 1000.0, 10.0);
    FILE *pf = fopen("./tests/test_eddy_mixing.txt", "w");

    for (i = 0; i < 100; ++i)
    {
        double val = 200.0 - 1.0 * (double)i;
        for (j = 0; j < 100; ++j)
            set_mat_ind(C, i, j, val);
    }

    mat* strm = allocate_matNxN(100);
    MODEL_SETTINGS.strm_density = 80.0;
    strmf_stratocumulus_sym(x, z, 1.0, strm);
    mat* w = gradient_mat(strm, x, 1);
    mat* u = gradient_mat(strm, z, 0); scl_mat(u, -1.0);
    MODEL_SETTINGS.strm_density = 270.0;
    free_mat(strm);


    unsigned long N;
    arakawa_2d* grid = make_arakawa_2d(C, x, z, u, w, NULL, 5.0);
    for (N = 0; N < 3600; ++N)
    {
        mpadvec2d(grid, 4, MPDATA_NONOSCIL | MPDATA_CYCLIC_X | MPDATA_NILL_Y | MPDATA_NFLUX_Y);
       // for (i = 0; i < 100; ++i) 
        //{
        //    grid->data[0][i] = 200.0; grid->data[99][i] = 100.0;
        //}
    }
    memcpy_arakawa_2d(C, grid);
    free_arakawa_2d(grid, 0);

    for (i = 0; i < 100; ++i)
    {
        for (j = 0; j < 100; ++j)
            fprintf(pf, "%.3e ", C->data[i * 100 + j]);
        fprintf(pf, "\r\n");
    }

    fclose(pf);
    free_mat(u); free_mat(w);
    free_vec(x); free_vec(z);
    free_mat(C);

}