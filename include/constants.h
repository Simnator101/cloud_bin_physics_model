#ifndef _H_CONSTANTS
#define _H_CONSTANTS

#include <math.h>
#define RAIR 287.058    // J / kg / K
#define RV   461.5      // J / kg / K
#define CPA  1005.0     //
#define CVA  718.0      // 
#define GVAL 9.81       // m / s^2
#define LV   2.5e6    // J / kg 
#define RHOW 997.0      // kg / m^3
#define TRBE 4.0        // m^2 / s^2 Turbulent Kinetic Energy
#define KIN_AIR 1.5e-5  // m^2 / s
#define KELVIN_C 1.2e-9 // m

#ifndef PI
#define PI   3.141592653589793
#endif

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

_CCODE_START
typedef enum __interp_type
{
    LINEAR,
    LOG
} interp_type;

/*Saturarion Calculations*/
double sat_press_vapour(const double T);
double sat_mixr_vapour(const double p, const double T);
double dewpoint_temp(const double p, const double q);

_CCODE_END

#endif