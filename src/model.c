#include <stdio.h>
#include <time.h>

#include "../include/vector.h"
#include "../include/pgrid.h"
#include "../include/environment.h"
#include "../include/mathfuncs.h"
#include "../include/extio.h"
#include "../include/mpdata.h"

#include "../include/tests.h"

int main(int argc, char **argv)
{
    // Tests
    //test_advection();
    //test_advection_2d();
    //test_droplet_advc();
    //test_droplet_coal();
    //test_eddy_mixing();
    //return 0;

    // Indices
    long long n, nf;
    register unsigned long i, j, k;
    unsigned long ind;
    unsigned long nc_i;
    unsigned long NZ, NX, NR;
    double t, fdt;
    //
    double Ncn; 
    double S;
    double Crate;
    // Flags
    int cyclic_diff_calc = 0;
    // Data Fields
    vec* x, *z, *r, *dr;
    mat* p0;
    mat* T, *Th, *Th_conv, *T_prev;
    mat* q, *l;
    mat* strmfnc, *rhow, *rhou, *rho, *rhovT_d;
    mat* Nc;
    pgrid* bins;
    mat* ccoll, *cck;
    vec* drdt;
    vec* Th0_bot, *q0_bot;
    vec* Th0_left, *q0_left;
    vec* Th0_top, *q0_top;
    vec* Th0_right, *q0_right;
    int* ima;
    // Progress bar
    unsigned pbi;
    float pdone = -1.0;
    char pbar[60]; memset(pbar, 0, sizeof(pbar));
    // Clock
    clock_t clck;
    clck = clock();

    // Load in options... 
    assert(argc > 1);
    assert(read_job_settings(argv[1]) == 0);
    fprint_opts(stdout, MODEL_SETTINGS);

    // If Cyclic in x set he horizontal diff calc flag
    cyclic_diff_calc = (MODEL_SETTINGS.zx_border_type & MPDATA_CYCLIC_X);

    // Initialise Spatial Scales
    x = range(MODEL_SETTINGS.xMin,
              MODEL_SETTINGS.xMax + MODEL_SETTINGS.dx,
              MODEL_SETTINGS.dx);
    z = range(MODEL_SETTINGS.zMin,
              MODEL_SETTINGS.zMax + MODEL_SETTINGS.dz,
              MODEL_SETTINGS.dz);
    switch (MODEL_SETTINGS.bgrid)
    {
    case LINEXP_RGRID:
        r = linexp_grid(MODEL_SETTINGS.nbins, MODEL_SETTINGS.dr_lin, MODEL_SETTINGS.dr_ex);
        break;
    case LINMASS_RGRID:
        r = linmass_grid(MODEL_SETTINGS.nbins, MODEL_SETTINGS.dr_lin, MODEL_SETTINGS.dr_ex);
        break;
    default:
        return -1;
    }
    dr = dr_bins(r);
    NZ = LEN(z); NX = LEN(x); NR = LEN(r);

    // Initialise Momentum Framework
    strmfnc = zero_mat(NZ, NX);
    rhow   = zero_mat(NZ, NX);
    rhou   = zero_mat(NZ, NX);

    // Initialise Hydrological Variables
    T    = allocate_matNxM(NZ, NX);
    Th   = allocate_matNxM(NZ, NX);
    q    = allocate_matNxM(NZ, NX);
    l    = zero_mat(NZ, NX);
    bins = zero_pgrid(NZ, NX, NR); 
    Nc = zero_mat(NZ, NX);

    // Load data from sounding to load in initial T
    sounding snd = {NULL, NULL, NULL, NULL};
    read_sounding_txt(&snd, MODEL_SETTINGS.snd_file_nm, MODEL_SETTINGS.snd_fmt);
    for (i = 0; i < LEN(z); ++i)
    {
        double _Thf = pow(snd.pe->data[0] / sample_snd_p(z->data[i], snd), RAIR / CPA);

        double _Tv  = sample_snd_Te(z->data[i], snd);
        double _qv  = sample_snd_qe(z->data[i], snd);
        //printf("Z %.2f m T %.2f K q %.2f g/kg\n", z->data[i], _Tv, _qv * 1e3);
        for (j = 0; j < LEN(x); ++j)
        {
            set_mat_ind(T, i, j, _Tv);
            set_mat_ind(Th, i, j, _Tv * _Thf);
            set_mat_ind(q, i, j, _qv);
        }
    }
    //mafprint_vec(stdout, snd.pe); fprint_vec(stdout, snd.ze);
    MODEL_SETTINGS.p0 = snd.pe->data[0];
    MODEL_SETTINGS.z0 = snd.ze->data[0];
    Th0_bot = slice_mat_row(Th, 0); 
    Th0_top = slice_mat_row(Th, NZ - 1); 
    Th0_left = slice_mat_col(Th, 0);
    Th0_right = slice_mat_col(Th, NX-1);

    q0_bot  = slice_mat_row(q, 0);
    q0_top  = slice_mat_row(q, NZ - 1);
    q0_left = slice_mat_col(q, 0);
    q0_right = slice_mat_col(q, NX-1);
    free_sounding(&snd);

    // Calculate Anelastic Background Profile
    // dummy_rho = one_vec(LEN(r));
    rho = ideal_gass_density(T, z, MODEL_SETTINGS.p0);
    p0 = ideal_gass_pressure(T, z, MODEL_SETTINGS.p0);
    Th_conv = potential_T_converter(rho, T);

    // TEST set rho to one
    // for (i = 0; i < LEN_MAT(rho); ++i) rho->data[i] = 1.0;

    // Initialise Droplet V therminal momentum
    rhovT_d = allocate_matNxM(NZ, LEN(r));
    for (i = 0; i < NZ; ++i) for (j = 0; j < LEN(r); ++j)
        set_mat_ind(rhovT_d, i, j, rho->data[i * NX] * droplet_vT(r->data[j]));
    drdt = droplet_growth_rate(r);
    
    // Calculate Collision Kernel
    ima = malloc(sizeof(int) * LEN(r) * LEN(r));
    ccoll = allocate_matNxN(LEN(r));
    courant_coal(r, ccoll, ima);
    cck = collision_kernel(r, KERNEL_LONG, MODEL_SETTINGS.dt);

    // Keep previous T
    T_prev = copy_mat(T);

    // Create Output file
    netCDF_obj out = init_netcdf_output(MODEL_SETTINGS.nc_output_nm, z, x, r);
    assert(out.state_var == 0);


    // Run Simulation
    for (nc_i = 0, n = 0; n < MODEL_SETTINGS.NT; ++n)
    {
        t = (double)n * MODEL_SETTINGS.dt;
        fdt = MODEL_SETTINGS.dt / (double)MODEL_SETTINGS.fNT;

        // Update Kinematic Framework
        switch (MODEL_SETTINGS.strm_type)
        {
        case STRM_SHALLOW_CUMULUS:
            strmf_shallow_cumulus(x, z, t, strmfnc);
            break;
        case STRM_SYM_EDDY:
            strmf_stratocumulus_sym(x, z, t, strmfnc);
            break;
        case STRM_ASYM_EDDY:
            strmf_stratocumulus_asym(x, z, t, strmfnc);
            break;
        }        
        free_mat(rhow); free_mat(rhou);
        rhow = gradient_mat(strmfnc, x, 1);
        rhou = gradient_mat(strmfnc, z, 0); scl_mat(rhou, -1.0);

        // Model Assumptions
        // No vertical Fluxes Through the upper and lower boundary
        memset(rhow->data + (NZ - 1)*NX, 0, sizeof(double)*NX);
        memset(rhow->data, 0, sizeof(double)*NX);
        // Also a no-slip condition on both vertical boundaries
        //memset(rhou->data + (NZ - 1)*NX, 0, sizeof(double) * NX);
        //memset(rhou->data, 0, sizeof(double) * NX);


        // With Eddies there is a large possibility of NaN development in the botom and top layer
        // To fix this, we clear these regions of ANY liquid
        //memset(q->data, 0, sizeof(double) * NX);
        //memset(q->data + (NZ - 1) * NX, 0, sizeof(double) * NX);
        //memset(bins->data, 0, NR * NX * sizeof(double));
        //memset(bins->data + (NZ-1) * NR * NX, 0, NR * NX * sizeof(double));
        
        // Write To netCDF    
        if ((n % MODEL_SETTINGS.nc_write_freq) == 0)
        {
            netcdf_add_t_entry(out, nc_i, t);
            write_netcdf_field(out, out.T_field, nc_i, T);
            write_netcdf_field(out, out.Th_field, nc_i, Th);
            write_netcdf_field(out, out.q_field, nc_i, q);
            write_netcdf_field(out, out.l_field, nc_i, l);
            write_netcdf_field(out, out.rhou_field, nc_i, rhou);
            write_netcdf_field(out, out.rhow_field, nc_i, rhow);
            write_netcdf_binvar(out, out.bins_field, nc_i, bins);
            ++nc_i;
        }


        // Fluxes due to SHF and LHF
        if (MODEL_SETTINGS.LHF != 0.0)
        {
            SHF_T_flux(T, rho, max_vec(z), MODEL_SETTINGS.dt);
            for(i = 0; i < LEN_MAT(T); ++i) Th->data[i] = Th_conv->data[i] * T->data[i];
            LHF_q_flux(q, rho, max_vec(z), MODEL_SETTINGS.dt);
        }

        // Advection of Temperature and Humidity
        arakawa_2d* Th_grid = make_arakawa_2d(Th, x, z, rhou, rhow, rho, MODEL_SETTINGS.dt);
        arakawa_2d* q_grid  = make_arakawa_2d(q, x, z, rhou, rhow, rho, MODEL_SETTINGS.dt);
        mpadvec2d(Th_grid, 2, MODEL_SETTINGS.zx_border_type);
        mpadvec2d(q_grid, 2, MODEL_SETTINGS.zx_border_type);
        memcpy_arakawa_2d(Th, Th_grid);
        memcpy_arakawa_2d(q, q_grid);
        free_arakawa_2d(Th_grid, 0); free_arakawa_2d(q_grid, 0); 
        

        // mpadvec(Th, cu, cw, rho, 2, BOUNDARY_CONTINUOUS | NON_OSCILLATORY);
        // mpadvec(q, cu, cw, rho, 2, BOUNDARY_CONTINUOUS | NON_OSCILLATORY);
        for (k = 0; k < NR; ++k)
        {
            mat* slc_b = slice_pgrid_x(bins, k);
            mat* wb = copy_mat(rhow);
            for (i = 0; i < NZ; ++i) for (j = 0; j < NX; ++j)
                wb->data[i*NX+j] = rhow->data[i*NX+j] - rhovT_d->data[i*NR+k];
            
            arakawa_2d* bin_grid = make_arakawa_2d(slc_b, x, z, rhou, wb, rho, MODEL_SETTINGS.dt);
            mpadvec2d(bin_grid, 2, MODEL_SETTINGS.zx_border_type);
            memcpy_arakawa_2d(slc_b, bin_grid);
            free_arakawa_2d(bin_grid, 0);
            set_pgrid_x(bins, slc_b, k);

            free_mat(slc_b);
            free_mat(wb);
        }
        

        // Horizontal Diffusion
        for (i = 0; i < NZ; ++i)
        {
            vec* Th_slc = slice_mat_row(Th, i);
            vec* rho_slc = slice_mat_row(rho, i);
            horizontal_diffusion(Th_slc, rho_slc, x, MODEL_SETTINGS.dt, cyclic_diff_calc);
            set_mat_row(Th, Th_slc, i);

            vec* q_slc = slice_mat_row(q, i);
            horizontal_diffusion(q_slc, rho_slc, x, MODEL_SETTINGS.dt, cyclic_diff_calc);
            set_mat_row(q, q_slc, i);

            for (k = 0; k < LEN(r); ++k)
            {
                vec* slc_l = slice_pgrid_zx(bins, i, k);
                horizontal_diffusion(slc_l, rho_slc, x, MODEL_SETTINGS.dt, cyclic_diff_calc);
                set_pgrid_zx(bins, slc_l, i, k);
                free_vec(slc_l);
            }

            free_vec(Th_slc); free_vec(q_slc); free_vec(rho_slc);
        }

        // Fix Edges
        //set_mat_row(q, q0_bot, 0); set_mat_row(q, q0_bot, 0);
        //set_mat_row(Th, Th0_bot, 0); set_mat_row(q, q0_bot, 0);
        //set_mat_row(Th, Th0_top, NZ - 1); set_mat_row(q, q0_bot, 0);;
        //set_mat_col(Th, Th0, NX-1); set_mat_col(q, q0, NX-1);
        //for (j = 0; j < NX; ++j) for (k = 0; k < LEN(r); ++k)
        //    bins->data[j * bins->NX + k] = 0.0;

        // Solve Temperature Equation from pot. T.
        memcpy(T->data, Th->data, sizeof(double) * LEN_MAT(T));
        div_mat(T, 1, Th_conv);

        
        // Solve Bin equations (2nd Part)
        // We already calculated advection + diffusion in z,x space
        // Next we solve the remaining right-hand equations, which are source terms
        // We then use smaller fractional time steps to do condensational advection
        //
        for (i = 0; i < NZ; ++i)
        {
            for (j = 0; j < NX; ++j)
            {
                ind = i * NX + j;
                vec* bin = slice_pgrid_zy(bins, i, j);
                // Activation
                Ncn = activated_CCN(q->data[ind], T->data[ind], p0->data[ind]);
                double dCCN = MAX(Ncn - Nc->data[ind], 0.0) / MODEL_SETTINGS.dt;
                bin->data[0] += dCCN / dr->data[0];
                if (Ncn > Nc->data[ind]) Nc->data[ind] = Ncn;

                // Coalescene/Collision
                bin = collision(bin, r, ima, ccoll, cck);    
                clip_vec(bin, 0.0, INFINITY);            
                
                // Condensation/Evaporation
                for (nf = 0; nf < MODEL_SETTINGS.fNT; ++nf)
                {
                    // Calculate Droplet Growth velocity (m/s)
                    // and its courant velocity
                    S = (q->data[ind] / sat_mixr_vapour(p0->data[ind], T->data[ind])) - 1.0;
                    vec* grwth = courant_bin(drdt, dr, S, fdt);

                    // Create Advection Object
                    vec* bn = advec_bin(bin, grwth);
                    clip_vec(bn, 0.0, INFINITY);

                    for (k = 0; k < NR; ++k)
                    {
                        grwth->data[k] = (bn->data[k] - bin->data[k]) / fdt * droplet_mass(r->data[k]);
                    }
                    Crate = trapz(grwth, r);
                    memcpy(bin->data, bn->data, LEN(bn) * sizeof(double));

                    // Apply Rates
                    T->data[ind] = T->data[ind] + LV / CPA * Crate * fdt;
                    q->data[ind] = q->data[ind] - Crate * fdt;
                    Th->data[ind] = T->data[ind] * Th_conv->data[ind];

                    // Free vec
                    free_vec(grwth);
                    free_vec(bn);

                }
                set_pgrid_zy(bins, bin, i, j);
                // Update Water Liquid Field
                for (k = 0; k < LEN(r); ++k)
                    bin->data[k] *= droplet_mass(r->data[k]);
                l->data[i * NX + j] = trapz(bin, r);
                free_vec(bin);
            }
        }

        // If mode is quasi-equilibrium we break if balance is reached
        // We Ignore this for now
        /*
        if (MODEL_SETTINGS.NT == INT64_MAX)
        {
            double MSE = 0.0;
            for (i = 0; i < LEN_MAT(T); ++i)
                MSE = MSE + (T_prev->data[i] - T->data[i]) * (T_prev->data[i] - T->data[i]);
            //if (MSE <= QUASI_EQUILIBRIUM_SENSITIVITY) MODEL_SETTINGS.NT = 0;
            //printf("MSE %.1e\r\n", MSE);
            if (MSE < MSE_prev) printf("\n%.4e\n", MSE);
            MSE_prev = MSE;
            
        }
        */
        
        // Progress bar only shown if run length is already known
        float nprog = (float)(n + 1) / (float)(MODEL_SETTINGS.NT);
        if ((nprog > pdone) && (MODEL_SETTINGS.NT != INT64_MAX))
        {
            pbar[0] = '\r'; 
            pbar[1] = '|'; 

            for (pbi = 2; pbi - 2 < (unsigned)(50.0f * nprog); ++pbi)
                pbar[pbi] = '#';
            for (; pbi - 2 < 50; ++pbi)
                pbar[pbi] = ' ';
            sprintf(pbar + pbi, "| %03u%%", (unsigned)(nprog * 100.0f));
            printf(pbar);
            fflush(stdout);

            pdone = nprog;
        }     

        // Save old matrix
        memcpy(T_prev->data, T->data, LEN_MAT(T) * sizeof(double));   
    }
    printf("\n\r");

    // Write final NC entry
    t = (double)MODEL_SETTINGS.NT * MODEL_SETTINGS.dt;
    netcdf_add_t_entry(out, nc_i, t);
    write_netcdf_field(out, out.T_field, nc_i, T);
    write_netcdf_field(out, out.Th_field, nc_i, Th);
    write_netcdf_field(out, out.q_field, nc_i, q);
    write_netcdf_field(out, out.l_field, nc_i, l);
    write_netcdf_field(out, out.rhow_field, nc_i, rhow);
    write_netcdf_binvar(out, out.bins_field, nc_i, bins);
    ++nc_i;
    netcdf_close(out);

    // Free Vectors
    free_vec(x); free_vec(z); free_vec(r); free_vec(dr);
    free_vec(drdt); 
    free_vec(Th0_bot); free_vec(q0_bot);
    free_vec(Th0_left); free_vec(q0_left);
    free_vec(Th0_top); free_vec(q0_top);
    free_vec(Th0_right); free_vec(q0_right);

    // Free Matrices
    free_mat(strmfnc); free_mat(rhow); free_mat(rhou); free_mat(rho);
    free_mat(p0);
    free_mat(rhovT_d);
    free_mat(T); free_mat(Th); free_mat(Th_conv); free_mat(T_prev);
    free_mat(q); free_mat(l);
    free_mat(Nc); 
    free_mat(ccoll); free_mat(cck);

    // Free Pgrids
    free_pgrid(bins);

    // Free Other
    free_model_settings();
    free(ima);

    // Clock
    clck = clock() - clck;

    // Print Resource sum
    printf("Simulation done --- Run Final Info:\n");
    printf("\tTotal Vector Memory Leaks: %lli\n", ACTIVE_VECTORS);
    printf("\tTotal Matrix Memory Leaks: %lli\n", ACTIVE_MATRICES);
    printf("\tTotal p-grid Memory Leaks: %lli\n", ACTIVE_PGRIDS);
    printf("\tElapsed Time: %.2f s\n", ((double)clck) / CLOCKS_PER_SEC);

    return 0;
}