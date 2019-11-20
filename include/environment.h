#ifndef _H_ENVIRONMENT
#define _H_ENVIRONMENT

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include <stdio.h>
#include "../include/mathfuncs.h"
#include "./constants.h"
#include "./vector.h"
#include "./matrix.h"

// Sounding Plain TXT format
#define SOUNDING_T      0x01U
#define SOUNDING_TD     0x02U
#define SOUNDING_P      0x04U
#define SOUNDING_Z      0x08U
#define SOUNDING_Q      0x10U
#define SOUNDING_TH     0x20U
#define SOUNDING_END    0x00U

#define KERNEL_LONG     0x01
#define KERNEL_HALL     0x02
#define KERNEL_GOLV     0x03

_CCODE_START

/*Conversion Calculations*/
vec* convert_z_to_p(vec* z, vec* T, double p0);
vec* convert_p_to_z(vec* p, vec* T, double z0);

/*Enviroment Calculations*/
mat* ideal_gass_pressure(mat* T0, vec* z0, const double p0);
mat* ideal_gass_density(mat* T0, vec* z0, const double p0);
mat* potential_T_converter(mat* rho, mat* T);

/*Sounding Module*/
typedef struct __sounding
{
    vec* pe;
    vec* ze;
    vec* Te;
    vec* Tde;
} sounding;

void free_sounding(sounding* snd);
int read_sounding_txt(sounding* snd, const char* file_nm, unsigned long long fmt);

#define sample_snd_Te(z, snd) interp(z, snd.ze, snd.Te, LINEAR)
#define sample_snd_Tde(z, snd) interp(z, snd.ze, snd.Tde, LINEAR)
#define sample_snd_p(z, snd) interp(z, snd.ze, snd.pe, LINEAR)
double sample_snd_qe(double z, sounding snd);

/*Kinematic Framework*/
void strmf_shallow_cumulus(vec* x, vec* y, double t, mat* strm);
void strmf_stratocumulus_sym(vec* x, vec* y, double t, mat* strm);
void strmf_stratocumulus_asym(vec* x, vec* y, double t, mat* strm);

/*Droplet Physics*/
#define droplet_mass(r) (4. / 3. * PI * (r) * (r) * (r) * RHOW)
double droplet_vT(double r);
double ventilation_coef(double r);
vec* droplet_growth_rate(vec* r);
void courant_coal(vec* r, mat* c, int* ima);
mat* collision_kernel(vec* r, const int type, const double dt);
vec* collision(vec* bin, vec* r, int* ima, mat* c, mat* ck);

/*Bin Definition and Conversion Formulas*/
vec* linexp_grid(const unsigned long N, const double linT, const double expT);
vec* linmass_grid(const unsigned long N, const double linT, const double massS);
mat* edges(const vec* r);
vec* dr_bins(const vec* r);

/*Subscale Turbulence*/
vec* horizontal_diffusion(vec* v, vec* rho, vec* x, const double dt, int cyclic_field);

/*Surface Parametrisation Fluxes*/
mat* SHF_T_flux(mat* T, const mat* rho, const double Hs, const double dt);
mat* LHF_q_flux(mat* q, const mat* rho, const double Hs, const double dt);

/*Droplet Advection Physics*/
vec* courant_bin(vec* drdt, vec* dr, const double S, const double dt);
vec* advec_bin(vec* b, vec* C_drdt);


_CCODE_END
#endif