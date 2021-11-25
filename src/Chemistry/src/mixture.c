/*******************************************************************************
 * @file mixture.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <math.h>
#include "chemistry_private.h"

/*******************************************************************************
 * @brief Calculate the mixtures heat capacity (const. pressure)
 * @param Y
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_cp(double *Y, double T, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * calc_species_cp_r(i, T, chemistry) * Rsp[i];

    return tmp;
}

/*******************************************************************************
 * @brief Calculate the mixtures heat capacity (const. volume)
 * @param cp
 * @param R
 * @return double
 ******************************************************************************/
double calc_mix_cv(double cp, double R)
{
    return cp - R;
}

/*******************************************************************************
 * @brief Calculate the mixtures free Gibb's energy
 * @param h
 * @param s
 * @param T
 * @return double
 ******************************************************************************/
double calc_mix_g(double h, double s, double T)
{
    return h - s * T;
}

/*******************************************************************************
 * @brief Calculate the mixtures enthalpy
 * @param Y
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_h(double *Y, double T, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * calc_species_h_rt(i, T, chemistry) * Rsp[i] * T;

    return tmp;
}

/*******************************************************************************
 * @brief Calculate the mixtures molar mass from mole fraction
 * @param X
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_molar_mass_from_X(double *X, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += X[i] * molar_mass[i];

    return tmp;
}

/*******************************************************************************
 * @brief Calculate the mixtures molar mass from mass fraction
 * @param Y
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_molar_mass_from_Y(double *Y, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] / molar_mass[i];

    return (1.0 / tmp);
}

/*******************************************************************************
 * @brief Calculate the mixtures gas constant
 * @param Y
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_R(double *Y, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * Rsp[i];

    return tmp;
}

/*******************************************************************************
 * @brief Calculate the mixtures entropy (with Dalton's fix)
 * @param X
 * @param Y
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_s(double *X, double *Y,
                  double p, double T, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
    {
        double tmp2 = calc_species_s_r(i, T, chemistry) - log(p / P0);
        if (X[i] > 0.0)
            tmp2 -= log(X[i]);
        tmp += Y[i] * tmp2 * Rsp[i];
    }

    return tmp;
}

/*******************************************************************************
 * @brief Calculate the mixtures internal energy
 * @param h
 * @param R
 * @param T
 * @return double
 ******************************************************************************/
double calc_mix_u(double h, double R, double T)
{
    return h - R * T;
}

/*******************************************************************************
 * @brief Convert the give mass fractions to concentrations
 * @param Y
 * @param rho
 * @param C
 * @param chemistry
 ******************************************************************************/
void conv_mass_to_conc(double *Y, double rho,
                       double *C, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        C[i] = rho * Y[i] / molar_mass[i];
    }
}

/*******************************************************************************
 * @brief Convert the give mass fractions to mole fractions
 * @param Y
 * @param mix_molar_mass
 * @param X
 * @param chemistry
 ******************************************************************************/
void conv_mass_to_mole(double *Y, double mix_molar_mass,
                       double *X, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        X[i] = mix_molar_mass * Y[i] / molar_mass[i];
    }
}

/*******************************************************************************
 * @brief Convert the give mole fractions to concentrations
 * @param X
 * @param rho
 * @param molar_mass
 * @param C
 * @param chemistry
 ******************************************************************************/
void conv_mole_to_conc(double *X, double rho, double molar_mass,
                       double *C, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
    {
        C[i] = rho * X[i] / molar_mass;
    }
}

/*******************************************************************************
 * @brief Convert the give mole fractions to concentrations
 * @param X
 * @param mix_molar_mass
 * @param Y
 * @param chemistry
 ******************************************************************************/
void conv_mole_to_mass(double *X, double mix_molar_mass,
                       double *Y, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        Y[i] = X[i] * molar_mass[i] / mix_molar_mass;
    }
}

/*******************************************************************************
 * @brief Correct the given mass/mole fractions
 * @param Y
 * @param Y_out
 * @param chemistry
 ******************************************************************************/
void correct_fraction(double *Y, double *Y_out, chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;

    double Y_sum = 0.0;
    for (int i = 0; i < n_specii - 1; ++i)
    {
        Y_out[i] = MAX(0.0, MIN(1.0, Y[i]));
        Y_sum += Y_out[i];
    }

    Y_out[n_specii - 1] = MAX(0.0, 1.0 - Y_sum);
    Y_sum += Y_out[n_specii - 1];

    for (int i = 0; i < n_specii - 1; ++i)
    {
        Y_out[i] = Y_out[i] / Y_sum;
    }
}
