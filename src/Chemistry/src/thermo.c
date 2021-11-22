/*******************************************************************************
 * @file thermo.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <math.h>
#include "chemistry_private.h"

const double s12 = 1. / 2.;   /** Calculation constant */
const double s13 = 1. / 3.;   /** Calculation constant */
const double s14 = 1. / 4.;   /** Calculation constant */
const double s15 = 1. / 4.;   /** Calculation constant */
const double s16 = 1. / 6.;   /** Calculation constant */
const double s112 = 1. / 12.; /** Calculation constant */
const double s120 = 1. / 20.; /** Calculation constant */

/*******************************************************************************
 * @brief Calculate dimensionless heat capacity
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_cp_r(int i, double T, chemistry_t *chemistry)
{
    double *a = (T < chemistry->specii->bounds[BOUNDS * i + 1])
                    ? &chemistry->specii->coeff_low[NASA * i]
                    : &chemistry->specii->coeff_high[NASA * i];

    return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + a[4] * T)));
}

/*******************************************************************************
 * @brief Calculate dimensionless enthalpy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_h_rt(int i, double T, chemistry_t *chemistry)
{
    double *a = (T < chemistry->specii->bounds[BOUNDS * i + 1])
                    ? &chemistry->specii->coeff_low[NASA * i]
                    : &chemistry->specii->coeff_high[NASA * i];

    double tmp = (s14 * a[3] + s15 * a[4] * T);
    tmp = s12 * a[1] + T * (s13 * a[2] + T * tmp);
    return a[0] + T * tmp + a[5] / T;
}

/*******************************************************************************
 * @brief Calculate dimensionless entropy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_s_r(int i, double T, chemistry_t *chemistry)
{
    double *a = (T < chemistry->specii->bounds[BOUNDS * i + 1])
                    ? &chemistry->specii->coeff_low[NASA * i]
                    : &chemistry->specii->coeff_high[NASA * i];

    double tmp = s13 * a[3] + s14 * a[4] * T;
    tmp = a[1] + T * (s12 * a[2] + T * tmp);
    return a[0] * log(T) + T * tmp + a[6];
}

/*******************************************************************************
 * @brief Calculate dimensionless free gibb's energy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_g_rt(int i, double T, chemistry_t *chemistry)
{
    double *a = (T < chemistry->specii->bounds[BOUNDS * i + 1])
                    ? &chemistry->specii->coeff_low[NASA * i]
                    : &chemistry->specii->coeff_high[NASA * i];

    double tmp = s112 * a[3] + s120 * a[4] * T;
    tmp = s12 * a[1] + T * (s16 * a[2] + T * tmp);
    return a[0] * (log(T) - 1.0) + T * tmp - a[5] / T + a[6];
}