/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "analyze_module.h"
#include "reactor/reactor_module.h"

double *residual = NULL;

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void analyze_define()
{
    REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    REGISTER_FINALIZE_ROUTINE(analyze_finalize);
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void analyze_finalize()
{
    DEALLOCATE(residual);
}

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void analyze_initialize()
{
    residual = ALLOCATE(sizeof(residual) * n_variables);
    set_value_n(0.0, n_variables, residual);
}

/*******************************************************************************
 * @brief Calculate the global residual
 * @param dt
 ******************************************************************************/
void calc_global_residual(double dt)
{
    for (int i = 0; i < n_variables; ++i)
        residual[i] = ABS(phi_dt[i]) * dt;
}