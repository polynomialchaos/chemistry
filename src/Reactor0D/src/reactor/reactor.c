/*******************************************************************************
 * @file xxx.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "reactor_module.h"




void_reactor_fp_t reactor_function_pointer = NULL;

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


void reactor_initialize();
void reactor_finalize();

void calc_reactor_isobar_adiabt(double t);
void calc_reactor_isobar_isotherm(double t);
void calc_reactor_isochor_adiabat(double t);
void calc_reactor_isochor_isotherm(double t);
double calc_temp_equation(double *dY_dt, double T_old);


void reactor_define()
{
    REGISTER_INITIALIZE_ROUTINE(reactor_initialize);
    REGISTER_FINALIZE_ROUTINE(reactor_finalize);

    string_t tmp = "untitled";
    SET_PARAMETER("Reactor/chem_file", StringParameter, &tmp, "The chemistry data file", NULL, 0);

    string_t tmp1_opt[] = {"ISOBAR-ADIABAT", "ISOBAR-ISOTHERM", "ISOCHOR-ADIABAT", "ISOCHOR-ISOTHERM"};
    int tmp1_opt_n = sizeof(tmp1_opt) / sizeof(string_t);
    string_t tmp1 = tmp1_opt[0];
    SET_PARAMETER("Reactor/type", StringParameter, &tmp1, "The reactor type", tmp1_opt, tmp1_opt_n);
}

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
        check_abort(0);
    }

    global_chemistry = read_chemistry_data(chem_file);
    global_state = allocate_state(global_chemistry);

    // read the input parameters
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
        check_abort(0);

    update_state_isobaric(global_state->p, global_state->rho, global_state->T, global_state->Y, global_state);

    // initialize thermochemical vectors
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

void calc_reactor_isobar_adiabt(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isobaric(global_state->p, global_state->rho, phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = -calc_temp_equation(&phi_dt[i_Y0], global_state->T) / global_state->cp;
}

void calc_reactor_isobar_isotherm(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isobaric(global_state->p, global_state->rho, phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = 0.0;
}

void calc_reactor_isochor_adiabat(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isochoric(global_state->p, global_state->rho, phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = -calc_temp_equation(&phi_dt[i_Y0], global_state->T) / global_state->cv;
}

void calc_reactor_isochor_isotherm(double t)
{
#ifdef DEBUG
    UNUSED(t);
#endif /* DEBUG */

    update_state_isochoric(global_state->p, global_state->rho, phi[0], &phi[i_Y0], global_state);
    calc_production_rate(global_state->C, global_state->T, global_chemistry);

    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    for (int i = 0; i < n_specii; ++i)
        phi_dt[i_Y0 + i] = specii->omega[i] / global_state->rho;

    phi_dt[0] = 0.0;
}

double calc_temp_equation(double *dY_dt, double T_old)
{
    specii_t *specii = global_chemistry->specii;
    int n_specii = specii->n_specii;

    double tmp = 0.0;
    for (int i = 0; i < n_specii; ++i)
        tmp += dY_dt[i] * calc_species_h_rt(i, T_old, global_chemistry) * specii->Rsp[i] * T_old;

    return tmp;
}
