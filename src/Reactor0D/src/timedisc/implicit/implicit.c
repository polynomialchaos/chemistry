
/*******************************************************************************
 * @file implicit.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <math.h>
#include "implicit_module.h"
#include "reactor/reactor_module.h"
#include "timedisc/timedisc_module.h"

int implicit_active = 0;
int n_iter_inner = 0;
int n_iter_lsoe = 0;
double **phi_old = NULL;

string_t implicit_scheme_name = NULL;
string_t method_name = NULL;
int max_iter_inner = 100;
string_t solver_name = NULL;
string_t jacobian_type_name = NULL;
double tolerance_lsoe = 1e-3;
int max_iter_lsoe = 100;
int max_krylov_dims = 15;
int max_krylov_restarts = 2;

math_int_ft matrix_vector = NULL;
int is_bicgstab = 1;
int n_work_size = 0;
double dt_loc = 0.0;
double tpdt_loc = 0.0;

enum
{
    NBDFStagesEuler = 1,
    NBDFStagesBDF2 = 2
};

double bdf_a_euler[NBDFStagesEuler] = {0.0};
double bdf_b_euler = 1.0;

double bdf_a_bdf2[NBDFStagesBDF2] = {-1. / 3., 1 / 3.};
double bdf_b_bdf2 = 2. / 3.;

int n_bdf_stages = 0;
double *bdf_a = NULL;
double bdf_b = 1.0;

int n_bdf_stages_loc = 0;
double *bdf_a_loc = NULL;
double bdf_b_loc = 1.0;

double *work = NULL;

double *Y_n = NULL;
double *f_Y_n = NULL;
double *dY_n = NULL;
double *dY_dt_n = NULL;
double *jac = NULL;

/*******************************************************************************
 * @brief Numerical jacobian routine
 * @param n_var
 ******************************************************************************/
void calc_jacobian_numerical(int n_var)
{
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
        double eps_fd = 0.0;
        eps_fd = Y_n[i_var] * Y_n[i_var];
        eps_fd = sqrt(eps_fd) * 1e-4;

        // positive + eps
        for (int j = 0; j < n_var; ++j)
            phi[j] = Y_n[j];

        phi[i_var] += 0.5 * eps_fd;
        phi[i_var] = MAX(phi_bounds[i_var * BOUNDDIM],
                         MIN(phi_bounds[i_var * BOUNDDIM + 1], phi[i_var]));

        reactor_function_pointer(tpdt_loc);

        for (int j = 0; j < n_var; ++j)
        {
            jac[i_var * n_var + j] = -phi_dt[j];
        }

        // negative + eps
        for (int j = 0; j < n_var; ++j)
            phi[j] = Y_n[j];

        phi[i_var] -= 0.5 * eps_fd;
        phi[i_var] = MAX(phi_bounds[i_var * BOUNDDIM],
                         MIN(phi_bounds[i_var * BOUNDDIM + 1], phi[i_var]));

        reactor_function_pointer(tpdt_loc);

        for (int j = 0; j < n_var; ++j)
        {
            jac[i_var * n_var + j] += phi_dt[j];
            jac[i_var * n_var + j] *= bdf_b_loc / (eps_fd + SMALL);
        }

        jac[i_var * n_var + i_var] += 1.0 / dt_loc;
    }
}

/*******************************************************************************
 * @brief Define implicit timedisc
 ******************************************************************************/
void implicit_define()
{
    REGISTER_INITIALIZE_ROUTINE(implicit_initialize);
    REGISTER_FINALIZE_ROUTINE(implicit_finalize);

    string_t tmp_opt[] = {"BDF-2", "Euler"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    string_t tmp2_opt[] = {"Newton"};
    int tmp2_opt_n = sizeof(tmp2_opt) / sizeof(string_t);
    string_t tmp2 = tmp2_opt[0];

    string_t tmp3_opt[] = {"Numerical"};
    int tmp3_opt_n = sizeof(tmp3_opt) / sizeof(string_t);
    string_t tmp3 = tmp3_opt[0];

    string_t tmp4_opt[] = {"BiCGStab", "GMRes"};
    int tmp4_opt_n = sizeof(tmp4_opt) / sizeof(string_t);
    string_t tmp4 = tmp4_opt[0];

    SET_PARAMETER("TimeDisc/Implicit/scheme", StringParameter, &tmp,
                  "The implicit timestep scheme", &tmp_opt, tmp_opt_n);
    SET_PARAMETER("TimeDisc/Implicit/method", StringParameter, &tmp2,
                  "The method to solve the non-linear system of equations",
                  &tmp2_opt, tmp2_opt_n);
    SET_PARAMETER("TimeDisc/Implicit/max_iter_inner", DigitParameter,
                  &max_iter_inner,
                  "The maximum number of inner iterations", NULL, 0);
    SET_PARAMETER("TimeDisc/Implicit/solver", StringParameter, &tmp4,
                  "The linear system of equations solver", &tmp4_opt,
                  tmp4_opt_n);
    SET_PARAMETER("TimeDisc/Implicit/jacobian_type", StringParameter, &tmp3,
                  "The type of jacobian generation", &tmp3_opt, tmp3_opt_n);
    SET_PARAMETER("TimeDisc/Implicit/tolerance_lsoe", NumberParameter,
                  &tolerance_lsoe,
                  "The linear solver tolerance", NULL, 0);
    SET_PARAMETER("TimeDisc/Implicit/max_iter_lsoe", DigitParameter,
                  &max_iter_lsoe,
                  "The linear solver maximum number of iterations", NULL, 0);
    SET_PARAMETER("TimeDisc/Implicit/max_krylov_dims", DigitParameter,
                  &max_krylov_dims,
                  "The maximum Krylov space dimension in GMRes solver", NULL, 0);
    SET_PARAMETER("TimeDisc/Implicit/max_krylov_restarts", DigitParameter,
                  &max_krylov_restarts,
                  "The maximum restarts performed in GMRes solver", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize implicit timedisc
 ******************************************************************************/
void implicit_finalize()
{
    DEALLOCATE(implicit_scheme_name);
    DEALLOCATE(method_name);
    DEALLOCATE(solver_name);
    DEALLOCATE(jacobian_type_name);

    for (int i = 0; i < n_bdf_stages; ++i)
        DEALLOCATE(phi_old[i]);
    DEALLOCATE(phi_old);

    DEALLOCATE(work);

    DEALLOCATE(Y_n);
    DEALLOCATE(f_Y_n);
    DEALLOCATE(dY_n);
    DEALLOCATE(dY_dt_n);
    DEALLOCATE(jac);
}

/*******************************************************************************
 * @brief Initialize implicit timedisc
 ******************************************************************************/
void implicit_initialize()
{
    if (implicit_active == 0)
        return;

    GET_PARAMETER("TimeDisc/Implicit/scheme", StringParameter,
                  &implicit_scheme_name);
    GET_PARAMETER("TimeDisc/Implicit/method", StringParameter,
                  &method_name);
    GET_PARAMETER("TimeDisc/Implicit/max_iter_inner", DigitParameter,
                  &max_iter_inner);
    GET_PARAMETER("TimeDisc/Implicit/solver", StringParameter,
                  &solver_name);
    GET_PARAMETER("TimeDisc/Implicit/jacobian_type", StringParameter,
                  &jacobian_type_name);
    GET_PARAMETER("TimeDisc/Implicit/tolerance_lsoe", NumberParameter,
                  &tolerance_lsoe);
    GET_PARAMETER("TimeDisc/Implicit/max_iter_lsoe", DigitParameter,
                  &max_iter_lsoe);
    GET_PARAMETER("TimeDisc/Implicit/max_krylov_dims", DigitParameter,
                  &max_krylov_dims);
    GET_PARAMETER("TimeDisc/Implicit/max_krylov_restarts", DigitParameter,
                  &max_krylov_restarts);

    if (is_equal(implicit_scheme_name, "Euler"))
    {
        n_bdf_stages = NBDFStagesEuler;
        bdf_a = bdf_a_euler;
        bdf_b = bdf_b_euler;
    }
    else if (is_equal(implicit_scheme_name, "BDF-2"))
    {
        n_bdf_stages = NBDFStagesBDF2;
        bdf_a = bdf_a_bdf2;
        bdf_b = bdf_b_bdf2;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    if (is_equal(method_name, "Newton"))
    {
        time_step_function_pointer = time_step_newton;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    if (is_equal(solver_name, "BiCGStab"))
    {
        is_bicgstab = 1;
        n_work_size = get_bicgstab_n_m_work_size(n_variables, 1);
        work = ALLOCATE(sizeof(double) * n_work_size);
    }
    else if (is_equal(solver_name, "GMRes"))
    {
        is_bicgstab = 0;
        n_work_size = get_gmres_n_m_work_size(n_variables, 1, max_krylov_dims);
        work = ALLOCATE(sizeof(double) * n_work_size);
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    if (is_equal(jacobian_type_name, "Numerical"))
    {
        matrix_vector = matrix_vector_numerical;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    phi_old = ALLOCATE(sizeof(double *) * n_bdf_stages);
    for (int i = 0; i < n_bdf_stages; ++i)
        phi_old[i] = ALLOCATE(sizeof(double) * n_variables);

    Y_n = ALLOCATE(sizeof(double) * n_variables);
    f_Y_n = ALLOCATE(sizeof(double) * n_variables);
    dY_n = ALLOCATE(sizeof(double) * n_variables);
    dY_dt_n = ALLOCATE(sizeof(double) * n_variables);
    jac = ALLOCATE(sizeof(double) * n_variables * n_variables);
}

/*******************************************************************************
 * @brief Implicit time discretizazion routine (Newton)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_newton(int iter, double t, double dt)
{
    // discretization parameters
    n_bdf_stages_loc = (iter < n_bdf_stages) ? NBDFStagesEuler : n_bdf_stages;
    bdf_a_loc = (iter < n_bdf_stages) ? bdf_a_euler : bdf_a;
    bdf_b_loc = (iter < n_bdf_stages) ? bdf_b_euler : bdf_b;

    // set local timestep
    dt_loc = dt;
    tpdt_loc = t + dt_loc;

    // store the old state before calculating the FVTimeDerivative
    for (int i = n_bdf_stages_loc - 1; i > 0; --i)
        copy_n(phi_old[i - 1], n_variables, phi_old[i]);

    for (int j = 0; j < n_variables; ++j)
        phi_old[0][j] = phi[j];

    // fill inital values for newton iteration
    copy_n(phi_old[0], n_variables, Y_n);
    copy_n(phi_dt, n_variables, dY_dt_n);

    // calculate the inital error for newton abort criterion
    for (int j = 0; j < n_variables; ++j)
    {
        int idx = j;
        f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc +
                     bdf_b_loc * dY_dt_n[idx];
    }

    for (int i_stage = 0; i_stage < n_bdf_stages_loc; ++i_stage)
        for (int j = 0; j < n_variables; ++j)
        {
            int idx = j;
            f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc * phi_old[i_stage][idx];
        }

    const double err_f_Y_0 = len_n(f_Y_n, n_variables);
    double err_f_Y_old = err_f_Y_0;

    for (n_iter_inner = 1; n_iter_inner <= max_iter_inner; ++n_iter_inner)
    {
        calc_jacobian_numerical(n_variables);

        /* Jac * dY = fY_n => dY ... Jacobian is determined
        via finite difference. fY_n = phi - dt * RHS */
        n_iter_lsoe = max_iter_lsoe;
        double residual_lsoe = tolerance_lsoe * err_f_Y_old;

        if (is_bicgstab)
        {
            solve_bicgstab_n_m(n_variables, 1, f_Y_n, dY_n,
                               work, matrix_vector,
                               &n_iter_lsoe, &residual_lsoe);
        }
        else
        {
            solve_gmres_n_m(n_variables, 1, f_Y_n, dY_n,
                            work, matrix_vector, &n_iter_lsoe,
                            &residual_lsoe, max_krylov_dims,
                            max_krylov_restarts);
        }

        // Y^(n+1) = Y^(n) + (Y^(n+1)-Y^(n))
        for (int j = 0; j < n_variables; ++j)
        {
            int idx = j;
            Y_n[idx] += dY_n[idx];
            phi[j] = Y_n[idx];
        }

        reactor_function_pointer(tpdt_loc);
        copy_n(phi_dt, n_variables, dY_dt_n);

        for (int j = 0; j < n_variables; ++j)
        {
            int idx = j;
            f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc +
                         bdf_b_loc * dY_dt_n[idx];
        }

        for (int i_stage = 0; i_stage < n_bdf_stages_loc; ++i_stage)
            for (int j = 0; j < n_variables; ++j)
            {
                int idx = j;
                f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc *
                              phi_old[i_stage][idx];
            }

        err_f_Y_old = len_n(f_Y_n, n_variables);
        if (err_f_Y_old < err_f_Y_0)
            break;
        if (is_transient == 0)
            break;
    }

    if (n_iter_inner >= max_iter_inner)
        CHECK_EXPRESSION(0);
}

/*******************************************************************************
 * @brief Matrix vector routine (called by solver)
 * @param x
 * @param b
 * @param n_var
 * @param m
 * @return int
 ******************************************************************************/
int matrix_vector_numerical(double *x, double *b, size_t n_var, size_t m)
{
#ifdef DEBUG
    UNUSED(m);
#endif /* DEBUG */

    for (size_t j = 0; j < n_var; ++j)
    {
        b[j] = jac[j] * x[0];
    }

    for (size_t i_var = 1; i_var < n_var; ++i_var)
    {
        for (size_t j = 0; j < n_var; ++j)
        {
            b[j] += jac[i_var * n_var + j] * x[i_var];
        }
    }

    return 0;
}