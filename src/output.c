#include "../include/extio.h"
#define ERR(X,O) {printf("NetCDF Error: %s\n", nc_strerror(X)); O.state_var = X; return O;}


netCDF_obj init_netcdf_output(const char* file_nm, vec* z, vec* x, vec* r)
{
    netCDF_obj obj;
    int ncret;

    memset(&obj, 0, sizeof(obj));
    if ((ncret = nc_create(file_nm, NC_WRITE, &obj.ncid)))
        ERR(ncret, obj);

    // Create Dimensions
    if ((ncret = nc_def_dim(obj.ncid, "time", NC_UNLIMITED, obj.dimids)))
        ERR(ncret, obj);
    if ((ncret = nc_def_dim(obj.ncid, "z", LEN(z), obj.dimids + 1)))
        ERR(ncret, obj);
    if ((ncret = nc_def_dim(obj.ncid, "x", LEN(x), obj.dimids + 2)))
        ERR(ncret, obj);
    if ((ncret = nc_def_dim(obj.ncid, "r", LEN(r), obj.dimids + 3)))
        ERR(ncret, obj);


    // Define Variables
    if ((ncret = nc_def_var(obj.ncid, "time", NC_FLOAT, 1, obj.dimids, obj.dimvars)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "z", NC_FLOAT, 1, obj.dimids + 1, obj.dimvars + 1)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "x", NC_FLOAT, 1, obj.dimids + 2, obj.dimvars + 2)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "r", NC_FLOAT, 1, obj.dimids + 3, obj.dimvars + 3)))
        ERR(ncerr, obj);

    // Assign Units
    if ((ncret = nc_put_att_text(obj.ncid, obj.dimvars[0], "units", 1, "s")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.dimvars[1], "units", 1, "m")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.dimvars[2], "units", 1, "m")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.dimvars[3], "units", 8, "10**-6 m")))
        ERR(ncerr, obj);

    //int surf_f_dims[2] = {obj.dimids[0], obj.dimids[2]};

    // Define Fields
    if ((ncret = nc_def_var(obj.ncid, "T", NC_FLOAT, 3, obj.dimids, &obj.T_field)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "Th", NC_FLOAT, 3, obj.dimids, &obj.Th_field)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "q", NC_FLOAT, 3, obj.dimids, &obj.q_field)))
        ERR(ncerr, obj);
    if ((ncret = nc_def_var(obj.ncid, "l", NC_FLOAT, 3, obj.dimids, &obj.l_field)))
        ERR(ncerr, obj);

    if ((ncret = nc_def_var(obj.ncid, "ld", NC_FLOAT, 4, obj.dimids, &obj.bins_field)))
        ERR(ncerr, obj);

    //if ((ncret = nc_def_var(obj.ncid, "precip", NC_FLOAT, 2, surf_f_dims, &obj.precip_amount)))
    //    ERR(ncerr, obj);


    // Optional Fields
    obj.rhow_field = 0;
    if (MODEL_SETTINGS.nc_additional_output & NC_OUTPUT_MOMENTUM)
    {
        if ((ncret = nc_def_var(obj.ncid, "rhow", NC_FLOAT, 3, obj.dimids, &obj.rhow_field)))
            ERR(ncerr, obj);           
        if ((ncret = nc_put_att_text(obj.ncid, obj.rhow_field, "units", 9, "kg/m**2/s")))
            ERR(ncerr, obj);
        if ((ncret = nc_put_att_text(obj.ncid, obj.rhow_field, "standard_name", 15, "momentum_flux_w")))
            ERR(ncerr, obj);
    }

    // Field Meta
    if ((ncret = nc_put_att_text(obj.ncid, obj.T_field, "units", 1, "K")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.T_field, "standard_name", 12, "temperature")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.Th_field, "units", 1, "K")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.Th_field, "standard_name", 21, "potential_temperature")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.q_field, "units", 5, "kg/kg")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.q_field, "standard_name", 19, "vapour_mixing_ratio")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.l_field, "units", 5, "kg/kg")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.l_field, "standard_name", 25, "water_liquid_mixing_ratio")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.bins_field, "units", 6, "1/m**4")))
        ERR(ncerr, obj);
    if ((ncret = nc_put_att_text(obj.ncid, obj.bins_field, "standard_name", 17, "water_bin_density")))
        ERR(ncerr, obj);

    //if ((ncret = nc_put_att_text(obj.ncid, obj.precip_amount, "units", 7, "kg/kg/m")))
    //    ERR(ncerr, obj);
    //if ((ncret = nc_put_att_text(obj.ncid, obj.precip_amount, "standard_name", 30, "total_precipitation_at_surface")))
    //    ERR(ncerr, obj);

    // End Def mode
    if ((ncret = nc_enddef(obj.ncid)))
        ERR(ncret, obj);

    // Write Coordinate system
    nc_put_var_double(obj.ncid, obj.dimvars[1], z->data);
    nc_put_var_double(obj.ncid, obj.dimvars[2], x->data);
    scl_vec(r, 1e6);    
    nc_put_var_double(obj.ncid, obj.dimvars[3], r->data);
    scl_vec(r, 1e-6);

    return obj;
}

void write_netcdf_field(netCDF_obj ncf, const int varid, const size_t ti, mat* m)
{
    if (varid != 0)
    {
        size_t start_p[3] = {ti, 0, 0,};
        size_t count_p[3] = {1, ROWS(m), COLS(m)};
        nc_put_vara_double(ncf.ncid, varid, start_p, count_p, m->data);
    }
}

void write_netcdf_binvar(netCDF_obj ncf, const int varid, const size_t ti, pgrid* p)
{
    if (varid != 0)
    {
        size_t start_p[4] = {ti, 0, 0, 0};
        size_t count_p[4] = {1, p->NZ, p->NY, p->NX};
        nc_put_vara_double(ncf.ncid, varid, start_p, count_p, p->data);
    }
}

void write_netcdf_sfluxvar(netCDF_obj ncf, const int varid, const size_t ti, vec* f)
{
    if (varid != 0)
    {
        size_t start_p[2] = {ti, 0};
        size_t count_p[2] = {1, f->N};
        nc_put_vara_double(ncf.ncid, varid, start_p, count_p, f->data);
    }
}