//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "analyze_module.h"
#include "reactor/reactor_module.h"




double *residual = NULL;


void analyze_initialize();
void analyze_finalize();


void analyze_define()
{
    register_initialize_routine(analyze_initialize);
    register_finalize_routine(analyze_finalize);
}

void analyze_initialize()
{
    residual = allocate(sizeof(residual) * n_variables);
    set_value_n(0.0, residual, n_variables);
}

void analyze_finalize()
{
    deallocate(residual);
}

void calc_global_residual(double dt)
{
    for (int i = 0; i < n_variables; ++i)
        residual[i] = u_abs(phi_dt[i]) * dt;
}