/*******************************************************************************
 * @file state.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "chemistry_private.h"

double calc_ig_p(double rho, double R, double T);
double calc_ig_rho(double p, double R, double T);
double calc_ig_T(double p, double rho, double R);

/*******************************************************************************
 * @brief Allocate the state structure
 * @param chemistry
 * @return state_t*
 ******************************************************************************/
state_t *allocate_state(chemistry_t *chemistry)
{
    state_t *state = ALLOCATE(sizeof(state_t));
    state->chemistry = chemistry;

    state->Y = ALLOCATE(sizeof(double) * chemistry->specii->n_specii);
    state->X = ALLOCATE(sizeof(double) * chemistry->specii->n_specii);
    /* 3rd body concentration */
    state->C = ALLOCATE(sizeof(double) * (chemistry->specii->n_specii + 1));

    set_value_n(0.0, chemistry->specii->n_specii, state->Y);
    set_value_n(0.0, chemistry->specii->n_specii, state->X);
    set_value_n(0.0, chemistry->specii->n_specii + 1, state->C);

    return state;
}

/*******************************************************************************
 * @brief Calculate the ideal gas pressure
 * @param rho
 * @param R
 * @param T
 * @return double
 ******************************************************************************/
double calc_ig_p(double rho, double R, double T)
{
    return rho * R * T;
}

/*******************************************************************************
 * @brief Calculate the ideal gas density
 * @param p
 * @param R
 * @param T
 * @return double
 ******************************************************************************/
double calc_ig_rho(double p, double R, double T)
{
    return p / (R * T);
}

/*******************************************************************************
 * @brief Calculate the ideal gas temperature
 * @param p
 * @param rho
 * @param R
 * @return double
 ******************************************************************************/
double calc_ig_T(double p, double rho, double R)
{
    return p / (rho * R);
}

/*******************************************************************************
 * @brief Deallocate the state structure
 * @param state
 ******************************************************************************/
void deallocate_state(state_t *state)
{
    if (state == NULL)
        return;

    DEALLOCATE(state->Y);
    DEALLOCATE(state->C);
    DEALLOCATE(state->X);
}

/*******************************************************************************
 * @brief Print the state structure
 * @param state
 ******************************************************************************/
void print_state(state_t *state)
{
    if (state == NULL)
        return;

    PRINTF("\n");
    printf_r_sep_title('-', "State");

    PRINTF("%s = %e\n", "p", state->p);
    PRINTF("%s = %e\n", "T", state->T);
    PRINTF("%s = %e\n", "rho", state->rho);

    for (int i = 0; i < state->chemistry->specii->n_specii; ++i)
    {
        if (state->Y[i] < YSMALL)
            continue;
        PRINTF("%s = %e (C=%e, X=%e)\n", state->chemistry->specii->symbol[i],
               state->Y[i], state->C[i], state->X[i]);
    }

    PRINTF("%s = %e\n", "R", state->R);
    PRINTF("%s = %e\n", "Molar mass", state->molar_mass);

    PRINTF("%s = %e\n", "cp", state->cp);
    PRINTF("%s = %e\n", "cv", state->cv);
    PRINTF("%s = %e\n", "h", state->h);
    PRINTF("%s = %e\n", "s", state->s);
    PRINTF("%s = %e\n", "u", state->u);
    PRINTF("%s = %e\n", "g", state->g);

    printf_r_sep('-');
    PRINTF("\n");
}

/*******************************************************************************
 * @brief Update the state structure (isochoric)
 * @param p
 * @param rho
 * @param T
 * @param Y
 * @param state
 ******************************************************************************/
void update_state_isochoric(double p, double rho, double T,
                            double *Y, state_t *state)
{
#ifdef DEBUG
    UNUSED(p);
#endif /* DEBUG */

    chemistry_t *chemistry = state->chemistry;

    state->rho = rho;
    state->T = T;

    correct_fraction(Y, state->Y, chemistry);

    state->R = calc_mix_R(state->Y, chemistry);
    state->p = calc_ig_p(state->rho, state->R, state->T);

    conv_mass_to_conc(state->Y, state->rho, state->C, chemistry);

    state->molar_mass = calc_mix_molar_mass_from_Y(state->Y, chemistry);

    conv_mass_to_mole(state->Y, state->molar_mass, state->X, chemistry);

    state->cp = calc_mix_cp(state->Y, state->T, chemistry);
    state->cv = calc_mix_cv(state->cp, state->R);
    state->h = calc_mix_h(state->Y, state->T, chemistry);
    state->s = calc_mix_s(state->X, state->Y, state->p, state->T, chemistry);
    state->u = calc_mix_u(state->h, state->R, state->T);
    state->g = calc_mix_g(state->h, state->s, state->T);
}

/*******************************************************************************
 * @brief Update the state structure (isobaric)
 * @param p
 * @param rho
 * @param T
 * @param Y
 * @param state
 ******************************************************************************/
void update_state_isobaric(double p, double rho, double T,
                           double *Y, state_t *state)
{
#ifdef DEBUG
    UNUSED(rho);
#endif /* DEBUG */

    chemistry_t *chemistry = state->chemistry;

    state->p = p;
    state->T = T;

    correct_fraction(Y, state->Y, chemistry);

    state->R = calc_mix_R(state->Y, chemistry);
    state->rho = calc_ig_rho(state->p, state->R, state->T);

    conv_mass_to_conc(state->Y, state->rho, state->C, chemistry);

    state->molar_mass = calc_mix_molar_mass_from_Y(state->Y, chemistry);

    conv_mass_to_mole(state->Y, state->molar_mass, state->X, chemistry);

    state->cp = calc_mix_cp(state->Y, state->T, chemistry);
    state->cv = calc_mix_cv(state->cp, state->R);
    state->h = calc_mix_h(state->Y, state->T, chemistry);
    state->s = calc_mix_s(state->X, state->Y, state->p, state->T, chemistry);
    state->u = calc_mix_u(state->h, state->R, state->T);
    state->g = calc_mix_g(state->h, state->s, state->T);
}