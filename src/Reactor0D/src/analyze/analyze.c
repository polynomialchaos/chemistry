/*******************************************************************************
 * @file analyze.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "analyze_module.h"
#include "reactor/reactor_module.h"

double *residual = NULL;

void analyze_initialize();
void analyze_finalize();

void analyze_define()
{
    REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    REGISTER_FINALIZE_ROUTINE(analyze_finalize);
}

void analyze_initialize()
{
    residual = ALLOCATE(sizeof(residual) * n_variables);
    set_value_n(0.0, n_variables, residual);
}

void analyze_finalize()
{
    DEALLOCATE(residual);
}

void calc_global_residual(double dt)
{
    for (int i = 0; i < n_variables; ++i)
        residual[i] = ABS(phi_dt[i]) * dt;
}