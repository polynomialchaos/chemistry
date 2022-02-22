/*******************************************************************************
 * @file reactor.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "reactor_module.h"

void_reactor_ft reactor_function_pointer = NULL;

string_t chem_file = NULL;
string_t reactor_type_name = NULL;

chemistry_t *global_chemistry = NULL;
state_t *global_state = NULL;

const int i_Y0 = 1;
int n_variables = 0;
string_t *variables = NULL;

double *phi = NULL;
double *phi_dt = NULL;
double *phi_bounds = NULL;

/*******************************************************************************
 * @brief Calculate the isobaric adiabatic 0D reactor
 * @param t
 ******************************************************************************/
void calc_reactor_isobar_adiabt(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isobaric(global_state->p, global_state->rho,
                          phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = -calc_dT_dt(&phi_dt[i_Y0], global_state->T) /
                global_state->cp;
}

/*******************************************************************************
 * @brief Calculate the isobaric isothermal 0D reactor
 * @param t
 ******************************************************************************/
void calc_reactor_isobar_isotherm(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isobaric(global_state->p, global_state->rho,
                          phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = 0.0;
}

/*******************************************************************************
 * @brief Calculate the isochoric adiabatic 0D reactor
 * @param t
 ******************************************************************************/
void calc_reactor_isochor_adiabat(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isochoric(global_state->p, global_state->rho,
                           phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = -calc_dT_dt(&phi_dt[i_Y0], global_state->T) /
                global_state->cv;
}

/*******************************************************************************
 * @brief Calculate the isochoric isothermal 0D reactor
 * @param t
 ******************************************************************************/
void calc_reactor_isochor_isotherm(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isochoric(global_state->p, global_state->rho,
                           phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = 0.0;
}

/*******************************************************************************
 * @brief Return the temperature equation source term
 * @param dY_dt
 * @param T_old
 * @return double
 ******************************************************************************/
double calc_dT_dt(double *dY_dt, double T_old)
{
    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += dY_dt[i] * calc_species_h_rt(i, T_old, global_chemistry) *
               specii->Rsp[i] * T_old;

    return tmp;
}

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void reactor_define()
{
    REGISTER_INITIALIZE_ROUTINE(reactor_initialize);
    REGISTER_FINALIZE_ROUTINE(reactor_finalize);

    string_t tmp = "untitled";
    SET_PARAMETER("Reactor/chem_file", StringParameter, &tmp,
                  "The chemistry data file", NULL, 0);

    string_t tmp1_opt[] = {"ISOBAR-ADIABAT", "ISOBAR-ISOTHERM",
                           "ISOCHOR-ADIABAT", "ISOCHOR-ISOTHERM"};
    int tmp1_opt_n = sizeof(tmp1_opt) / sizeof(string_t);
    string_t tmp1 = tmp1_opt[0];
    SET_PARAMETER("Reactor/type", StringParameter, &tmp1,
                  "The reactor type", tmp1_opt, tmp1_opt_n);

    double p = 1e5;
    SET_PARAMETER("Reactor/p", NumberParameter, &p, NULL, NULL, 0);

    double T = 1200.0;
    SET_PARAMETER("Reactor/T", NumberParameter, &T, NULL, NULL, 0);

    double Y_CH4 = 0.05;
    SET_PARAMETER("Reactor/Y/CH4", NumberParameter, &Y_CH4, NULL, NULL, 0);

    double Y_O2 = 0.1;
    SET_PARAMETER("Reactor/Y/O2", NumberParameter, &Y_O2, NULL, NULL, 0);

    double Y_N2 = 0.85;
    SET_PARAMETER("Reactor/Y/N2", NumberParameter, &Y_N2, NULL, NULL, 0);
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void reactor_finalize()
{
    DEALLOCATE(chem_file);
    DEALLOCATE(reactor_type_name);

    deallocate_chemistry(global_chemistry);
    DEALLOCATE(global_chemistry);

    deallocate_state(global_state);
    DEALLOCATE(global_state);

    for (int i = 0; i < n_variables; ++i)
        DEALLOCATE(variables[i]);
    DEALLOCATE(variables);

    DEALLOCATE(phi);
    DEALLOCATE(phi_dt);
    DEALLOCATE(phi_bounds);
}

/*******************************************************************************
 * @brief Initialize reactor
 ******************************************************************************/
void reactor_initialize()
{
    GET_PARAMETER("Reactor/chem_file", StringParameter, &chem_file);
    GET_PARAMETER("Reactor/type", StringParameter, &reactor_type_name);

    if (is_equal(reactor_type_name, "ISOBAR-ADIABAT"))
    {
        reactor_function_pointer = calc_reactor_isobar_adiabt;
    }
    else if (is_equal(reactor_type_name, "ISOBAR-ISOTHERM"))
    {
        reactor_function_pointer = calc_reactor_isobar_isotherm;
    }
    else if (is_equal(reactor_type_name, "ISOCHOR-ADIABAT"))
    {
        reactor_function_pointer = calc_reactor_isochor_adiabat;
    }
    else if (is_equal(reactor_type_name, "ISOCHOR-ISOTHERM"))
    {
        reactor_function_pointer = calc_reactor_isochor_isotherm;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    global_chemistry = read_chemistry_data(chem_file);
    global_state = allocate_state(global_chemistry);

    /* read the input parameters */
    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    GET_PARAMETER("Reactor/p", NumberParameter, &global_state->p);
    GET_PARAMETER("Reactor/T", NumberParameter, &global_state->T);

    for (int i = 0; i < n_specii; ++i)
    {
        string_t tmp = allocate_strcat("Reactor/Y/", specii->symbol[i]);
        if (PARAMETER_EXISTS(tmp))
            GET_PARAMETER(tmp, NumberParameter, &global_state->Y[i]);
        DEALLOCATE(tmp);
    }

    double Y_sum = sum_n(global_state->Y, n_specii);
    if (ABS(Y_sum - 1.0) > YONE)
        CHECK_EXPRESSION(0);

    update_state_isobaric(global_state->p, global_state->rho,
                          global_state->T, global_state->Y, global_state);

    /* initialize thermochemical vectors */
    n_variables = i_Y0 + specii->n_specii;
    variables = ALLOCATE(sizeof(string_t) * n_variables);

    variables[0] = allocate_strcpy("T");
    for (int i = 0; i < specii->n_specii; ++i)
        variables[i_Y0 + i] = allocate_strcpy(specii->symbol[i]);

    phi = ALLOCATE(sizeof(double) * n_variables);
    phi_dt = ALLOCATE(sizeof(double) * n_variables);
    phi_bounds = ALLOCATE(sizeof(double) * BOUNDDIM * n_variables);

    phi[0] = global_state->T;
    copy_n(global_state->Y, n_specii, &phi[i_Y0]);

    phi_bounds[0] = 300.0;
    phi_bounds[1] = 4500.0;
    for (int i = 0; i < n_specii; ++i)
    {
        phi_bounds[(i_Y0 + i) * BOUNDDIM] = 0.0;
        phi_bounds[(i_Y0 + i) * BOUNDDIM + 1] = 1.0;
    }

    print_state(global_state);
}
