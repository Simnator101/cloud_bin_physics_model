        /*
        // Solve Advection Equations Z-component
        for (i = 0; i < NX; ++i)
        {
            vec* slc_z = slice_mat_col(rhow, i);
            vec* cz = courant_bott(slc_z, z, MODEL_SETTINGS.dt);
            vec* rho_slc = slice_mat_col(rho, i);
            //vec* rho_slc = one_vec(NX);

            vec* slc_Th = slice_mat_col(Th, i);
            adv2p(slc_Th, cz, rho_slc, BOUNDARY_NILL);
            set_mat_col(Th, slc_Th, i);
            free_vec(slc_Th);

            vec* slc_q = slice_mat_col(q, i);
            adv2p(slc_q, cz, rho_slc, BOUNDARY_NILL);
            set_mat_col(q, slc_q, i);
            free_vec(slc_q);
            
            vec* qr_bflux = zero_vec(LEN(r));
            for (k = 0; k < LEN(r); ++k)
            {
                free_vec(cz);
                vec* slc_w = slice_mat_col(rhovT_d, k);
                sub_vec(slc_w, 1, slc_z); scl_vec(slc_w, -1.0);

                cz = courant_bott(slc_w, z, MODEL_SETTINGS.dt);
                vec* slc_l = slice_pgrid_yx(bins, i, k);
                slc_l = adv2p(slc_l, cz, rho_slc, BOUNDARY_NILL | BOUNDARY_KEEP_OUT_FLUXES);

                // Split off boundary fluxes
                qr_bflux->data[k] = slc_l->data[0] * droplet_mass(r->data[k]);
                pop_back(slc_l); pop_front(slc_l);

                set_pgrid_yx(bins, slc_l, i, k);
                free_vec(slc_l); free_vec(slc_w); 
            }
            
            // Calc Rain Density at bottom
            rain_am->data[i] += trapz(qr_bflux, r);
            free_vec(qr_bflux);

            free_vec(slc_z);
            free_vec(cz);    
            free_vec(rho_slc);      
        }

        // Solve Advection Equations + Dispersion Equations X-component
        for (j = 0; j < NZ; ++j)
        {
            vec* slc_x = slice_mat_row(rhou, j);
            vec* cx = courant_bott(slc_x, x, MODEL_SETTINGS.dt);
            vec* rho_slc = slice_mat_row(rho, j);

            vec* slc_Th = slice_mat_row(Th, j);
            adv2p(slc_Th, cx, rho_slc, BOUNDARY_CONTINUOUS);
            horizontal_diffusion(slc_Th, rho_slc, x, MODEL_SETTINGS.dt);
            set_mat_row(Th, slc_Th, j);
            free_vec(slc_Th);

            vec* slc_q = slice_mat_row(q, j);
            adv2p(slc_q, cx, rho_slc, BOUNDARY_CONTINUOUS);
            horizontal_diffusion(slc_q, rho_slc, x, MODEL_SETTINGS.dt);
            set_mat_row(q, slc_q, j);
            free_vec(slc_q);

            for (k = 0; k < LEN(r); ++k)          
            {
                vec* slc_l = slice_pgrid_zx(bins, j, k);
                adv2p(slc_l, cx, rho_slc, BOUNDARY_NILL);
                horizontal_diffusion(slc_l, rho_slc, x, MODEL_SETTINGS.dt);
                set_pgrid_zx(bins, slc_l, j, k);
                free_vec(slc_l);
            }

            free_vec(slc_x);
            free_vec(cx);
            free_vec(rho_slc);
        }
        */
