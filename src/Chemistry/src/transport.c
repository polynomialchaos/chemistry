/*******************************************************************************
 * @file transport.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <math.h>
#include "chemistry_private.h"

const double s16 = 1. / 6.;   /** Calculation constant */
const double s32 = 3. / 2.;   /** Calculation constant */
const double s38 = 3. / 8.;   /** Calculation constant */
const double s316 = 3. / 16.; /** Calculation constant */
const double s52 = 5. / 2.;   /** Calculation constant */
const double s53 = 5. / 3.;   /** Calculation constant */
const double s516 = 5. / 16.; /** Calculation constant */

/** Collision integral (11) constants */
const double omega_aj[] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};

/** Collision integral (22) constants */
const double omega_bj[] = {1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0};

/*******************************************************************************
 * @brief Calculate collision integral (11)
 * @param T_r
 * @param delta_k
 * @return double
 ******************************************************************************/
double calc_coll_int_11(double T_r, double delta_k)
{
    double tmp = (exp(omega_aj[4] / T_r) - exp(-omega_aj[5] / T_r));
    tmp = 1. + (tmp * delta_k * delta_k) / (2. + 2.5 * delta_k);
    return (omega_aj[0] * pow(T_r, -omega_aj[1]) +
            pow(T_r + omega_aj[2], -omega_aj[3])) *
           tmp;
}

/*******************************************************************************
 * @brief Calculate collision integral (22)
 * @param T_r
 * @param delta_k
 * @return double
 ******************************************************************************/
double calc_coll_int_22(double T_r, double delta_k)
{
    double tmp = (exp(omega_bj[4] / T_r) - exp(-omega_bj[5] / T_r));
    tmp = 1. + (tmp * delta_k * delta_k) / (2. + 2.5 * delta_k);
    return (omega_bj[0] * pow(T_r, -omega_bj[1]) +
            pow(T_r + omega_bj[2], -omega_bj[3])) *
           tmp;
}

/*******************************************************************************
 * @brief Calculate interaction parameter collision diameter
 * @param i
 * @param j
 * @param xi
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_inter_param_col_lj(int i, int j, double xi, chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;

    return 0.5 * (specii->col_lj[i] + specii->col_lj[j]) * pow(xi, -s16);
}

/*******************************************************************************
 * @brief Calculate interaction parameter molecule weight
 * @param i
 * @param j
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_inter_param_molecule_weight(int i, int j, chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;

    return (specii->molecule_weight[i] * specii->molecule_weight[j]) /
           (specii->molecule_weight[i] + specii->molecule_weight[j]);
}

/*******************************************************************************
 * @brief Calculate interaction parameter Lennard Jones potential
 * @param i
 * @param j
 * @param xi
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_inter_param_pot_lj(int i, int j, double xi, chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;

    return sqrt(specii->pot_lj[i] * specii->pot_lj[j]) * xi * xi;
}

/*******************************************************************************
 * @brief Calculate interaction parameter (xi)
 * @param i
 * @param j
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_inter_param_xi(int i, int j, chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;

    if ((specii->dip_mo[i] <= 0.0) && (specii->dip_mo[j] > 0.0))
    {
        double alpha_r = specii->pol[i] / (specii->col_lj[i] *
                                           specii->col_lj[i] *
                                           specii->col_lj[i]);
        double mu_r = specii->dip_mo[j] /
                      sqrt(specii->pot_lj[j] * KB * specii->col_lj[j] *
                           specii->col_lj[j] * specii->col_lj[j]);

        return 1. + 0.25 * alpha_r * mu_r *
                        sqrt(specii->pot_lj[j] / specii->pot_lj[i]);
    }
    else if ((specii->dip_mo[i] > 0.0) && (specii->dip_mo[j] <= 0.0))
    {
        double alpha_r = specii->pol[j] /
                         (specii->col_lj[j] * specii->col_lj[j] *
                          specii->col_lj[j]);
        double mu_r = specii->dip_mo[i] /
                      sqrt(specii->pot_lj[i] * KB * specii->col_lj[i] *
                           specii->col_lj[i] * specii->col_lj[i]);

        return 1. + 0.25 * alpha_r * mu_r *
                        sqrt(specii->pot_lj[i] / specii->pot_lj[j]);
    }

    return 1.0;
}

/*******************************************************************************
 * @brief Calculate rotational relaxation
 * @param T
 * @param pot_lj
 * @return double
 ******************************************************************************/
double calc_rot_rel(double T, double pot_lj)
{
    const double f1 = pow(MCPI, 1.5) / 2.;
    const double f2 = MCPI * MCPI / 4. + 2.;
    const double f3 = pow(MCPI, 1.5);

    double ratio = pot_lj / T;
    return 1. + f1 * sqrt(ratio) + f2 * ratio + f3 * pow(ratio, 1.5);
}

/*******************************************************************************
 * @brief Calculate species self diffusion coefficient
 * @param i
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_dii(int i, double p, double T,
                        chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;
    double col_lj_sq = specii->col_lj[i] * specii->col_lj[i];

    double T_r = T / specii->pot_lj[i];
    double delta_k = 0.5 * specii->dip_mo[i] * specii->dip_mo[i] /
                     (specii->pot_lj[i] * KB * col_lj_sq * specii->col_lj[i]);

    double tmp = p * MCPI * col_lj_sq * calc_coll_int_11(T_r, delta_k);
    return s38 * sqrt(MCPI * pow(KB * T, 3) / specii->molecule_weight[i]) / tmp;
}

/*******************************************************************************
 * @brief Calculate species binary diffusion coefficient
 * @param i
 * @param j
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_dij(int i, int j, double p, double T,
                        chemistry_t *chemistry)
{
    double xi = calc_inter_param_xi(i, j, chemistry);
    double pot_lj = calc_inter_param_pot_lj(i, j, xi, chemistry);
    double col_lj = calc_inter_param_col_lj(i, j, xi, chemistry);

    specii_t *specii = chemistry->specii;
    double col_lj_sq = col_lj * col_lj;

    double T_r = T / pot_lj;
    double delta_k = 0.5 * specii->dip_mo[i] * specii->dip_mo[j] /
                     (pot_lj * KB * col_lj_sq * col_lj);

    double tmp = p * MCPI * col_lj_sq * calc_coll_int_11(T_r, delta_k);
    double M = calc_inter_param_molecule_weight(i, j, chemistry);
    return s316 * sqrt(2. * MCPI * pow(KB * T, 3) / M) / tmp;
}

/*******************************************************************************
 * @brief Calculate species heat conductivity
 * @param i
 * @param p
 * @param T
 * @param mu
 * @param Dii
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_lambda(int i, double p, double T, double mu, double Dii,
                           chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;

    double cvT = 0.0;
    double cvR = 0.0;
    double cvV = 0.0;

    if (specii->geom[i] == 0) /* single atom */
    {
        cvT = s32 * RM;
        cvR = 0.0;
        cvV = 0.0;
    }
    else if (specii->geom[i] == 1) /* linear molecule */
    {
        cvT = s32 * RM;
        cvR = RM;
        cvV = (calc_species_cp_r(i, T, chemistry) - 1.) * RM -
              cvT - cvR; /* includes CV = CP - R */
    }
    else /* non-linear molcule */
    {
        cvT = s32 * RM;
        cvR = s32 * RM;
        cvV = (calc_species_cp_r(i, T, chemistry) - 1.) * RM -
              cvT - cvR; /* includes CV = CP - R */
    }

    double rho = p / (specii->Rsp[i] * T);
    double Z = specii->rot_rel[i] * calc_rot_rel(i, 298.0) / calc_rot_rel(i, T);
    double A = s52 - rho * Dii / mu;
    double B = Z + 2. / MCPI * (s53 * cvR / RM + rho * Dii / mu);

    double fT = s52 * (1 - 2. / MCPI * cvR / cvT * A / B);
    double fR = rho * Dii / mu * (1 + 2. / MCPI * A / B);
    double fV = rho * Dii / mu;

    return mu / specii->molar_mass[i] * (fT * cvT + fR * cvR + fV * cvV);
}

/*******************************************************************************
 * @brief Calculate species viscosity
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_mu(int i, double T, chemistry_t *chemistry)
{
    specii_t *specii = chemistry->specii;
    double col_lj_sq = specii->col_lj[i] * specii->col_lj[i];

    double T_r = T / specii->pot_lj[i];
    double delta_k = 0.5 * specii->dip_mo[i] * specii->dip_mo[i] /
                     (specii->pot_lj[i] * KB * col_lj_sq * specii->col_lj[i]);

    double tmp = MCPI * col_lj_sq * calc_coll_int_22(T_r, delta_k);
    return s516 * sqrt(MCPI * specii->molecule_weight[i] * KB * T) / tmp;
}