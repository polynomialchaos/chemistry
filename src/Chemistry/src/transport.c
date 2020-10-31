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
const double s16        = 1./6.;
const double s32        = 3./2.;
const double s38        = 3./8.;
const double s316       = 3./16.;
const double s52        = 5./2.;
const double s53        = 5./3.;
const double s516       = 5./16.;

const double omega_aj[] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};
const double omega_bj[] = {1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0};

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
double calc_rot_rel( double T, double pot_lj );
double calc_inter_param_molecule_weight( int i, int j, Chemistry_t *chemistry );
double calc_inter_param_pot_lj( int i, int j, double xi, Chemistry_t *chemistry );
double calc_inter_param_col_lj( int i, int j, double xi, Chemistry_t *chemistry );
double calc_inter_param_xi( int i, int j, Chemistry_t *chemistry );
double calc_coll_int_11( double T_r, double delta_k );
double calc_coll_int_22( double T_r, double delta_k );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
double calc_sp_mu( int i, double T, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;
    double col_lj_sq = specii->col_lj[i] * specii->col_lj[i];

    double T_r     = T / specii->pot_lj[i];
    double delta_k = 0.5 * specii->dip_mo[i] * specii->dip_mo[i] /
        (specii->pot_lj[i] * KB * col_lj_sq * specii->col_lj[i]);

    double tmp = PI * col_lj_sq * calc_coll_int_22( T_r, delta_k );
    return s516 * sqrt( PI * specii->molecule_weight[i] * KB * T ) / tmp;
}

double calc_sp_lambda( int i, double p, double T, double mu, double Dii, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;

    double cvT = 0.0;
    double cvR = 0.0;
    double cvV = 0.0;

    if (specii->geom[i] == 0) // single atom
    {
        cvT = s32 * RM;
        cvR = 0.0;
        cvV = 0.0;
    }
    else if (specii->geom[i] == 1) // linear molecule
    {
        cvT = s32 * RM;
        cvR = RM;
        cvV = (calc_sp_cp_r( i, T, chemistry ) - 1.) * RM - cvT - cvR; // includes CV = CP - R
    }
    else // non-linear molcule
    {
        cvT = s32 * RM;
        cvR = s32 * RM;
        cvV = (calc_sp_cp_r( i, T, chemistry ) - 1.) * RM - cvT - cvR; // includes CV = CP - R
    }

    double rho = p / (specii->Rsp[i] * T);
    double Z   = specii->rot_rel[i] * calc_rot_rel( i, 298.0 ) / calc_rot_rel( i, T );
    double A   = s52 - rho * Dii / mu;
    double B   = Z + 2. / PI * ( s53 * cvR / RM + rho * Dii / mu );

    double fT  = s52 * ( 1 - 2. / PI * cvR / cvT * A / B );
    double fR  = rho * Dii / mu * ( 1 + 2. / PI * A / B );
    double fV  = rho * Dii / mu;

    return mu / specii->molar_mass[i] * ( fT * cvT + fR * cvR + fV * cvV );
}

double calc_sp_dii( int i, double p, double T, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;
    double col_lj_sq = specii->col_lj[i] * specii->col_lj[i];

    double T_r     = T / specii->pot_lj[i];
    double delta_k = 0.5 * specii->dip_mo[i] * specii->dip_mo[i] /
        (specii->pot_lj[i] * KB * col_lj_sq * specii->col_lj[i]);

    double tmp = p * PI * col_lj_sq * calc_coll_int_11( T_r, delta_k );
    return s38 * sqrt( PI * pow( KB * T, 3 ) / specii->molecule_weight[i] ) / tmp;
}

double calc_sp_dij( int i, int j, double p, double T, Chemistry_t *chemistry )
{
    double xi       = calc_inter_param_xi( i, j, chemistry );
    double pot_lj   = calc_inter_param_pot_lj( i, j, xi, chemistry );
    double col_lj   = calc_inter_param_col_lj( i, j, xi, chemistry );

    Specii_t *specii = chemistry->specii;
    double col_lj_sq = col_lj * col_lj;

    double T_r      = T / pot_lj;
    double delta_k  = 0.5 * specii->dip_mo[i] * specii->dip_mo[j] / (pot_lj * KB * col_lj_sq * col_lj);

    double tmp = p * PI * col_lj_sq * calc_coll_int_11( T_r, delta_k );
    return s316 * sqrt( 2. * PI * pow( KB * T, 3 ) / calc_inter_param_molecule_weight( i, j, chemistry ) ) / tmp;
}

double calc_rot_rel( double T, double pot_lj )
{
    const double f1 = pow( PI, 1.5 ) / 2.;
    const double f2 = PI * PI / 4. + 2.;
    const double f3 = pow( PI, 1.5 );

    double ratio = pot_lj / T;
    return 1. + f1 * sqrt( ratio ) + f2 * ratio + f3 * pow( ratio, 1.5 );
}

double calc_inter_param_molecule_weight( int i, int j, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;

    return (specii->molecule_weight[i] * specii->molecule_weight[j]) /
        (specii->molecule_weight[i] + specii->molecule_weight[j]);
}

double calc_inter_param_pot_lj( int i, int j, double xi, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;

    return sqrt( specii->pot_lj[i] * specii->pot_lj[j] ) * xi * xi;
}

double calc_inter_param_col_lj( int i, int j, double xi, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;

    return 0.5 * (specii->col_lj[i] + specii->col_lj[j]) * pow( xi, -s16 );
}

double calc_inter_param_xi( int i, int j, Chemistry_t *chemistry )
{
    Specii_t *specii = chemistry->specii;

    if ((specii->dip_mo[i] <= 0.0) && (specii->dip_mo[j] > 0.0))
    {
        double alpha_r = specii->pol[i] / (specii->col_lj[i] * specii->col_lj[i] * specii->col_lj[i]);
        double mu_r = specii->dip_mo[j] /
            sqrt( specii->pot_lj[j] * KB * specii->col_lj[j] * specii->col_lj[j] * specii->col_lj[j] );

        return 1. + 0.25 * alpha_r * mu_r * sqrt( specii->pot_lj[j] / specii->pot_lj[i] );
    }
    else if ((specii->dip_mo[i] > 0.0) && (specii->dip_mo[j] <= 0.0))
    {
        double alpha_r = specii->pol[j] / (specii->col_lj[j] * specii->col_lj[j] * specii->col_lj[j]);
        double mu_r = specii->dip_mo[i] /
            sqrt( specii->pot_lj[i] * KB * specii->col_lj[i] * specii->col_lj[i] * specii->col_lj[i] );

        return 1. + 0.25 * alpha_r * mu_r * sqrt( specii->pot_lj[i] / specii->pot_lj[j] );
    }

    return 1.0;
}

double calc_coll_int_11( double T_r, double delta_k )
{
    double tmp = 1. + ((exp( omega_aj[4] / T_r ) - exp( -omega_aj[5] / T_r )) * delta_k * delta_k) / (2. + 2.5 * delta_k);
    return (omega_aj[0] * pow( T_r, -omega_aj[1] ) + pow( T_r + omega_aj[2], -omega_aj[3] )) * tmp;
}

double calc_coll_int_22( double T_r, double delta_k )
{
    double tmp = 1. + ((exp( omega_bj[4] / T_r ) - exp( -omega_bj[5] / T_r )) * delta_k * delta_k) / (2. + 2.5 * delta_k);
    return (omega_bj[0] * pow( T_r, -omega_bj[1] ) + pow( T_r + omega_bj[2], -omega_bj[3] )) * tmp;
}