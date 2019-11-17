#ifndef _H_EXTIO
#define _H_EXTIO

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <netcdf.h>

#include "./vector.h"
#include "./matrix.h"
#include "./pgrid.h"

#define OPTIONS_PRINT_UNKNOWN

#define LINEXP_RGRID 0
#define LINMASS_RGRID 1

#define NC_OUTPUT_MOMENTUM_U     0x01U
#define NC_OUTPUT_MOMENTUM_W     0x02U
#define NC_OUTPUT_CLOUD_LIQUID   0x04U

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define STR_FIND_FIRST 1
#define STR_FIND_LAST  2

_CCODE_START
/*General Tools*/
long file_size(FILE* pf);

/*Text Reader Tools*/ 
long str_find_char(const char* str, const char delim, int mode);
char* str_trim(char* str, int left, int right);
int str_eval_func(char* str, int(*func)(int));

typedef char** strbuffer;
unsigned long strbuffer_len(strbuffer buffer);
#define allocate_strbuffer(LEN) calloc((LEN) + 1, sizeof(char*))
void free_strbuffer(strbuffer strb);
strbuffer str_split(const char* str, const char delim);

/*Model Run Settings*/
typedef struct __type_model_settings
{
    // Time Scales
    unsigned long NT;
    double dt;
    unsigned int fNT;

    // Spatial Scales
    double zMin, zMax;
    double xMin, xMax;
    double dz, dx;

    // Bin Settings
    unsigned long nbins;
    int bgrid;
    double dr_lin, dr_ex;

    // Initial values (Used for conversions)
    double p0, z0;

    // Surface Fluxes
    int use_fluxes; double LHF;

    // Stream Function Settings
    double strm_density, Zclb, Ztop, Xwidth;

    // CCN Settings
    double CCN_C0, CCN_k;

    // Initial Sounding Data to set hydrological vars to
    const char* snd_file_nm;
    unsigned long long snd_fmt;

    // NetCDF Output
    const char* nc_output_nm;
    unsigned short nc_write_freq;
    unsigned long nc_flags;

    // Internal String Storage
    strbuffer __str_int;


} __model_settings;
extern __model_settings MODEL_SETTINGS;

char* model_settings_add_str(char* str);
#define free_model_settings() free_strbuffer(MODEL_SETTINGS.__str_int)
int read_job_settings(const char *file_name);
void fprint_opts(FILE* pf, __model_settings sets);

/*NetCDF output*/
typedef struct __netcdf_out
{
    int ncid;
    int dimids[4];
    int dimvars[4];

    int T_field, Th_field, q_field, l_field;
    int bins_field;

    int rhow_field, rhou_field;

    int state_var;
} netCDF_obj;

netCDF_obj init_netcdf_output(const char* file_nm, vec* z, vec* x, vec* r);
void write_netcdf_field(netCDF_obj ncf, const int varid, const size_t ti, mat* m);
void write_netcdf_binvar(netCDF_obj ncf, const int varid, const size_t ti, pgrid* p);
void write_netcdf_sfluxvar(netCDF_obj ncf, const int varid, const size_t ti, vec* f);
#define netcdf_add_t_entry(nobj, ind, t) nc_put_var1_double(nobj.ncid, nobj.dimvars[0], &ind, &t)
#define netcdf_close(OBJ) nc_close(OBJ.ncid)

_CCODE_END

#endif