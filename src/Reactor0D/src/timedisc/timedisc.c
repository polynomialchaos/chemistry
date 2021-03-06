/*******************************************************************************
 * @file timedisc.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "timedisc_module.h"
#include "explicit/explicit_module.h"
#include "implicit/implicit_module.h"
#include "restart/restart_module.h"
#include "analyze/analyze_module.h"
#include "output/output_module.h"
#include "reactor/reactor_module.h"

void_timestep_ft time_step_function_pointer = NULL;
double_calc_timestep_ft calc_time_step_function_pointer = NULL;

int is_viscous_dt = 0;

string_t time_step_name = NULL;
int max_iter = 1000;
int is_transient = 1;
double abort_residual = 1e-10;
double dt = 1e-6;
double t_start = 0.0;
double t_end = 0.2;

/*******************************************************************************
 * @brief Print the residual line
 * @param iter
 * @param t
 * @param dt
 * @param do_output
 ******************************************************************************/
void print_residual(int iter, double t, double dt, int do_output)
{
    char output_str = (do_output == 1) ? '*' : ' ';
    char viscous_str = (is_viscous_dt == 1) ? 'T' : 'F';

    if (explicit_active)
    {
        BM_PRINT("%09d %12.5e %12.5e %c %c:",
                 iter, t, dt, viscous_str, output_str);
    }
    else
    {
        BM_PRINT("%09d %12.5e %12.5e %c %c %6d %6d:",
                 iter, t, dt, viscous_str, output_str,
                 n_iter_inner, n_iter_lsoe);
    }

    BM_PRINT(" %12.5e", residual[0]);
    BM_PRINT(" %12.5e",
             sum_n(&residual[i_Y0], n_variables - 1) / (n_variables - 1));
    BM_PRINT(" %12.5e", max_n(&residual[i_Y0], n_variables - 1));

    BM_PRINT(" %12.5e", phi[0]);

    BM_PRINT("\n");
}

/*******************************************************************************
 * @brief Print the residual line header
 ******************************************************************************/
void print_residual_header()
{
    if (explicit_active)
    {
        BM_PRINT("%9s %12s %12s %1s %1s:",
                 "iter", "time", "dt", "V", "O");
    }
    else
    {
        BM_PRINT("%9s %12s %12s %1s %1s %6s %6s:",
                 "iter", "time", "dt", "V", "O", "inner", "lsoe");
    }

    BM_PRINT(" %12s", "T");
    BM_PRINT(" %12s", "Y_mean");
    BM_PRINT(" %12s", "Y_max");

    BM_PRINT(" %12s", "T");

    BM_PRINT("\n");
}

/*******************************************************************************
 * @brief Time discretizazion routine
 ******************************************************************************/
void timedisc()
{
    /* check for errors */
    check_abort(0);

    int do_output = 0;
    int do_finalize = 0;

    int iter = 0;
    if (use_restart == 1)
        iter = iter_restart;

    double t = t_start;
    if (use_restart == 1)
        t = t_restart;

    print_residual_header();
    if (iter == 0)
        reactor_function_pointer(t);
    if ((do_output_data == 1) && (use_restart == 0))
        write_output(iter, t);

    while (1)
    {
        /* check for errors */
        check_abort(0);

        /* calculate time step dt */
        if (t + dt > t_end)
        {
            dt = t_end - t;
            do_output = 1;
            do_finalize = 1;
        }

        /* the timestep to be called */
        time_step_function_pointer(iter, t, dt);
        calc_global_residual(dt);

        /* check for NAN and INF */
        if (is_nan_n(residual, n_variables) ||
            is_inf_n(residual, n_variables))
            BM_CHECK_EXPRESSION(0);

        t = t + dt;
        iter = iter + 1;

        /* steady-state simulation */
        if ((is_transient == 0) &&
            (min_n(residual, n_variables) < abort_residual))
        {
            t_end = t;
            do_output = 1;
            do_finalize = 1;
        }

        if ((i_output_data > 0) && (iter % i_output_data == 0))
        {
            do_output = 1;
        }

        /* maximum iteration number reacher */
        if (iter >= iter_restart + max_iter)
        {
            t_end = t;
            do_output = 1;
            do_finalize = 1;
        }

        print_residual(iter, t, dt, do_output);

        /* output */
        if ((do_output_data == 1) && (do_output == 1))
        {
            write_output(iter, t);
            do_output = 0;
        }

        /* stop if required */
        if (do_finalize)
            break;
    }
}

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(timedisc_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(timedisc_finalize);

    string_t tmp_opt[] = {"Explicit", "Implicit"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    BM_SET_PARAMETER("TimeDisc/time_step", StringParameter, &tmp,
                     "The timestep mehtod", &tmp_opt, tmp_opt_n);
    BM_SET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter,
                     "The maximum number of iterations", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient,
                     "The flag wheter to be transient or steady-state",
                     NULL, 0);
    BM_SET_PARAMETER("TimeDisc/abort_residual",
                     NumberParameter, &abort_residual,
                     "The abort residual", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/dt", NumberParameter, &dt,
                     "The timestep", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start,
                     "The start time", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end,
                     "The end time", NULL, 0);

    explicit_define();
    implicit_define();
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void timedisc_finalize()
{
    time_step_function_pointer = NULL;
    calc_time_step_function_pointer = NULL;

    BM_DEALLOCATE(time_step_name);
}

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize()
{
    BM_GET_PARAMETER("TimeDisc/time_step", StringParameter, &time_step_name);
    BM_GET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter);
    BM_GET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient);
    BM_GET_PARAMETER("TimeDisc/abort_residual",
                     NumberParameter, &abort_residual);
    BM_GET_PARAMETER("TimeDisc/dt", NumberParameter, &dt);
    BM_GET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start);
    BM_GET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end);

    if (is_transient == 0)
        t_end = BC_DMAX;

    if (is_equal(time_step_name, "Explicit"))
    {
        explicit_active = 1;
    }
    else if (is_equal(time_step_name, "Implicit"))
    {
        implicit_active = 1;
    }
    else
    {
        BM_CHECK_EXPRESSION(0);
    }
}
