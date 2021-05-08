//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <math.h>
#include "chemistry_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void correct_fraction(double *Y, double *Y_out, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;

    double Y_sum = 0.0;
    for (int i = 0; i < n_specii - 1; ++i)
    {
        Y_out[i] = u_max(0.0, u_min(1.0, Y[i]));
        Y_sum += Y_out[i];
    }

    Y_out[n_specii - 1] = u_max(0.0, 1.0 - Y_sum);
    Y_sum += Y_out[n_specii - 1];

    for (int i = 0; i < n_specii - 1; ++i)
    {
        Y_out[i] = Y_out[i] / Y_sum;
    }

    // if (u_abs( Y_sum - 1.0 ) > YONE)
    //     check_error( 0 );
}

double calc_mix_R(double *Y, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * Rsp[i];

    return tmp;
}

double calc_mix_molar_mass_from_X(double *X, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += X[i] * molar_mass[i];

    return tmp;
}

double calc_mix_molar_mass_from_Y(double *Y, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] / molar_mass[i];

    return (1.0 / tmp);
}

double calc_mix_cp(double *Y, double T, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * calc_sp_cp_r(i, T, chemistry) * Rsp[i];

    return tmp;
}

double calc_mix_cv(double cp, double R)
{
    return cp - R;
}

double calc_mix_h(double *Y, double T, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += Y[i] * calc_sp_h_rt(i, T, chemistry) * Rsp[i] * T;

    return tmp;
}

double calc_mix_s(double *X, double *Y, double p, double T, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *Rsp = chemistry->specii->Rsp;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
    {
        double tmp2 = calc_sp_s_r(i, T, chemistry) - log(p / P0);
        if (X[i] > 0.0)
            tmp2 -= log(X[i]);
        tmp += Y[i] * tmp2 * Rsp[i];
    }

    return tmp;
}

double calc_mix_u(double h, double R, double T)
{
    return h - R * T;
}

double calc_mix_g(double h, double s, double T)
{
    return h - s * T;
}

void conv_mass_to_mole(double *Y, double mix_molar_mass, double *X, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        X[i] = mix_molar_mass * Y[i] / molar_mass[i];
    }
}

void conv_mole_to_mass(double *X, double mix_molar_mass, double *Y, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        Y[i] = X[i] * molar_mass[i] / mix_molar_mass;
    }
}

void conv_mass_to_conc(double *Y, double rho, double *C, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;
    double *molar_mass = chemistry->specii->molar_mass;

    for (int i = 0; i < n_specii; ++i)
    {
        C[i] = rho * Y[i] / molar_mass[i];
    }
}

void conv_mole_to_conc(double *X, double rho, double molar_mass, double *C, Chemistry_t *chemistry)
{
    int n_specii = chemistry->specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
    {
        C[i] = rho * X[i] / molar_mass;
    }
}
