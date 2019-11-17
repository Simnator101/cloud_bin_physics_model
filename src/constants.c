#include "../include/constants.h"

double dewpoint_temp(const double p, const double q)
{
    const double T0 = 273.16;
    const double es0 = 610.78;
    const double pwf = RV * T0 / LV;

    double es = q * p / RAIR * RV;
    double expt = log(pow(es / es0, pwf));
    return T0 / (1. - expt + 1e-30);
}

double sat_press_vapour(const double T)
{
    const double T0 = 273.16;
    const double es0 = 610.78;
    return es0 * exp(LV / RV / T0 * (1. - T0 / T));
}

double sat_mixr_vapour(const double p, const double T)
{
    double es = sat_press_vapour(T);
    return es / p * RAIR / RV;
}