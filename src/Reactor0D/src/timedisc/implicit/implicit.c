
//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <math.h>
#include "implicit_module.h"
#include "reactor/reactor_module.h"
#include "timedisc/timedisc_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int implicit_active = 0;
int n_iter_inner    = 0;
int n_iter_lsoe     = 0;
double **phi_old    = NULL;

string_t implicit_scheme_name   = NULL;
string_t method_name            = NULL;
int max_iter_inner              = 100;
string_t solver_name            = NULL;
string_t jacobian_type_name     = NULL;
double tolerance_lsoe           = 1e-3;
int max_iter_lsoe               = 100;
int max_krylov_dims             = 15;
int max_krylov_restarts         = 2;

int_matvec_n_m_fp_t matrix_vector   = NULL;
int is_bicgstab                     = 1;
int n_work_size                     = 0;
double dt_loc                       = 0.0;
double tpdt_loc                     = 0.0;

enum {
    NBDFStagesEuler = 1,
    NBDFStagesBDF2  = 2
};

double bdf_a_euler[NBDFStagesEuler] = {0.0};
double bdf_b_euler                  = 1.0;

double bdf_a_bdf2[NBDFStagesBDF2]   = {-1./3.,1/3.};
double bdf_b_bdf2                   = 2./3.;

int n_bdf_stages    = 0;
double *bdf_a       = NULL;
double bdf_b        = 1.0;

int n_bdf_stages_loc    = 0;
double *bdf_a_loc       = NULL;
double bdf_b_loc        = 1.0;

double *work    = NULL;

double *Y_n     = NULL;
double *f_Y_n   = NULL;
double *dY_n    = NULL;
double *dY_dt_n = NULL;
double *jac     = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void implicit_initialize();
void implicit_finalize();

void calc_jacobian_numerical( int n_var );
int matrix_vector_numerical( double *x, double *b, int n, int m );
void time_step_newton( int iter, double t, double dt );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void implicit_define()
{
    register_initialize_routine( implicit_initialize );
    register_finalize_routine( implicit_finalize );

    string_t tmp_opt[] = {"BDF-2", "Euler"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    string_t tmp2_opt[] = {"Newton"};
    int tmp2_opt_n = sizeof( tmp2_opt ) / sizeof( string_t );
    string_t tmp2 = tmp2_opt[0];

    string_t tmp3_opt[] = {"Numerical"};
    int tmp3_opt_n = sizeof( tmp3_opt ) / sizeof( string_t );
    string_t tmp3 = tmp3_opt[0];

    string_t tmp4_opt[] = {"BiCGStab", "GMRes"};
    int tmp4_opt_n = sizeof( tmp4_opt ) / sizeof( string_t );
    string_t tmp4 = tmp4_opt[0];

    set_parameter( "TimeDisc/Implicit/scheme", ParameterString, &tmp,
        "The implicit timestep scheme", &tmp_opt, tmp_opt_n );
    set_parameter( "TimeDisc/Implicit/method", ParameterString, &tmp2,
        "The method to solve the non-linear system of equations", &tmp2_opt, tmp2_opt_n );
    set_parameter( "TimeDisc/Implicit/max_iter_inner", ParameterDigit, &max_iter_inner,
        "The maximum number of inner iterations", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/solver", ParameterString, &tmp4,
        "The linear system of equations solver", &tmp4_opt, tmp4_opt_n );
    set_parameter( "TimeDisc/Implicit/jacobian_type", ParameterString, &tmp3,
        "The type of jacobian generation", &tmp3_opt, tmp3_opt_n );
    set_parameter( "TimeDisc/Implicit/tolerance_lsoe", ParameterNumber, &tolerance_lsoe,
        "The linear solver tolerance", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_iter_lsoe", ParameterDigit, &max_iter_lsoe,
        "The linear solver maximum number of iterations", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_krylov_dims", ParameterDigit, &max_krylov_dims,
        "The maximum Krylov space dimension in GMRes solver", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_krylov_restarts", ParameterDigit, &max_krylov_restarts,
        "The maximum restarts performed in GMRes solver", NULL, 0 );
}

void implicit_initialize()
{
    if (implicit_active == 0) return;

    get_parameter( "TimeDisc/Implicit/scheme", ParameterString, &implicit_scheme_name );
    get_parameter( "TimeDisc/Implicit/method", ParameterString, &method_name );
    get_parameter( "TimeDisc/Implicit/max_iter_inner", ParameterDigit, &max_iter_inner );
    get_parameter( "TimeDisc/Implicit/solver", ParameterString, &solver_name );
    get_parameter( "TimeDisc/Implicit/jacobian_type", ParameterString, &jacobian_type_name );
    get_parameter( "TimeDisc/Implicit/tolerance_lsoe", ParameterNumber, &tolerance_lsoe );
    get_parameter( "TimeDisc/Implicit/max_iter_lsoe", ParameterDigit, &max_iter_lsoe );
    get_parameter( "TimeDisc/Implicit/max_krylov_dims", ParameterDigit, &max_krylov_dims );
    get_parameter( "TimeDisc/Implicit/max_krylov_restarts", ParameterDigit, &max_krylov_restarts );

    if (is_equal( implicit_scheme_name, "Euler" ))
    {
        n_bdf_stages    = NBDFStagesEuler;
        bdf_a           = bdf_a_euler;
        bdf_b           = bdf_b_euler;
    }
    else if (is_equal( implicit_scheme_name, "BDF-2" ))
    {
        n_bdf_stages    = NBDFStagesBDF2;
        bdf_a           = bdf_a_bdf2;
        bdf_b           = bdf_b_bdf2;
    }
    else
    {
        check_error( 0 );
    }

    if (is_equal( method_name, "Newton" ))
    {
        time_step_function_pointer = time_step_newton;
    }
    else
    {
        check_error( 0 );
    }

    if (is_equal( solver_name, "BiCGStab" ))
    {
        is_bicgstab = 1;
        n_work_size = get_bicgstab_n_m_work_size( n_variables, 1 );
        work = allocate( sizeof( double ) * n_work_size );
    }
    else if (is_equal( solver_name, "GMRes" ))
    {
        is_bicgstab = 0;
        n_work_size = get_gmres_n_m_work_size( n_variables, 1, max_krylov_dims );
        work = allocate( sizeof( double ) * n_work_size );
    }
    else
    {
        check_error( 0 );
    }

    if (is_equal( jacobian_type_name, "Numerical" ))
    {
        matrix_vector = matrix_vector_numerical;
    }
    else
    {
        check_error( 0 );
    }

    phi_old = allocate( sizeof( double* ) * n_bdf_stages );
    for ( int i = 0; i < n_bdf_stages; i++ )
        phi_old[i] = allocate( sizeof( double ) * n_variables );

    Y_n     = allocate( sizeof( double ) * n_variables );
    f_Y_n   = allocate( sizeof( double ) * n_variables );
    dY_n    = allocate( sizeof( double ) * n_variables );
    dY_dt_n = allocate( sizeof( double ) * n_variables );
    jac     = allocate( sizeof( double ) * n_variables * n_variables );
}

void implicit_finalize()
{
    deallocate( implicit_scheme_name );
    deallocate( method_name);
    deallocate( solver_name );
    deallocate( jacobian_type_name );

    for ( int i = 0; i < n_bdf_stages; i++ )
        deallocate( phi_old[i] );
    deallocate( phi_old );

    deallocate( work );

    deallocate( Y_n );
    deallocate( f_Y_n );
    deallocate( dY_n );
    deallocate( dY_dt_n );
    deallocate( jac );
}

void time_step_newton( int iter, double t, double dt )
{
    // discretization parameters
    n_bdf_stages_loc    = (iter < n_bdf_stages) ? NBDFStagesEuler : n_bdf_stages;
    bdf_a_loc           = (iter < n_bdf_stages) ? bdf_a_euler : bdf_a;
    bdf_b_loc           = (iter < n_bdf_stages) ? bdf_b_euler : bdf_b;

    // set local timestep
    dt_loc      = dt;
    tpdt_loc    = t + dt_loc;

    // store the old state before calculating the FVTimeDerivative
    for ( int i = n_bdf_stages_loc - 1; i > 0; i-- )
        copy_n( phi_old[i-1], phi_old[i], n_variables );

    for ( int j = 0; j < n_variables; j++ )
        phi_old[0][j] = phi[j];

    // fill inital values for newton iteration
    copy_n( phi_old[0], Y_n, n_variables );
    copy_n( phi_dt, dY_dt_n, n_variables );

    // calculate the inital error for newton abort criterion
    for ( int j = 0; j < n_variables; j++ )
    {
        int idx = j;
        f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc + bdf_b_loc * dY_dt_n[idx];
    }

    for ( int i_stage = 0; i_stage < n_bdf_stages_loc; i_stage++ )
        for ( int j = 0; j < n_variables; j++ )
        {
            int idx = j;
            f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc * phi_old[i_stage][idx];
        }

    const double err_f_Y_0 = len_n( f_Y_n, n_variables );
    double err_f_Y_old = err_f_Y_0;

    for ( n_iter_inner = 1; n_iter_inner <= max_iter_inner; n_iter_inner++ )
    {
        calc_jacobian_numerical( n_variables );

        // Jac * dY = fY_n => dY ... Jacobian is determined via finite difference. fY_n = phi - dt * RHS
        n_iter_lsoe = max_iter_lsoe;
        double residual_lsoe = tolerance_lsoe * err_f_Y_old;

        if (is_bicgstab)
        {
            BiCGStab_n_m( n_variables, 1, f_Y_n, dY_n,
                work, matrix_vector, &n_iter_lsoe, &residual_lsoe );
        }
        else
        {
            GMRes_n_m( n_variables, 1, f_Y_n, dY_n,
                work, matrix_vector, &n_iter_lsoe, &residual_lsoe, max_krylov_dims, max_krylov_restarts );
        }

        // Y^(n+1) = Y^(n) + (Y^(n+1)-Y^(n))
        for ( int j = 0; j < n_variables; j++ )
        {
            int idx = j;
            Y_n[idx] += dY_n[idx];
            phi[j] = Y_n[idx];
        }

        reactor_function_pointer( tpdt_loc );
        copy_n( phi_dt, dY_dt_n, n_variables );

        for ( int j = 0; j < n_variables; j++ )
        {
            int idx = j;
            f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc + bdf_b_loc * dY_dt_n[idx];
        }

        for ( int i_stage = 0; i_stage < n_bdf_stages_loc; i_stage++ )
            for ( int j = 0; j < n_variables; j++ )
            {
                int idx = j;
                f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc * phi_old[i_stage][idx];
            }

        err_f_Y_old = len_n( f_Y_n, n_variables );
        if (err_f_Y_old < err_f_Y_0) break;
        if (is_transient == 0) break;
    }

    if (n_iter_inner >= max_iter_inner)
        check_error( 0 );
}

void calc_jacobian_numerical( int n_var )
{
    for ( int i_var = 0; i_var < n_var; i_var++ )
    {
        double eps_fd = 0.0;
        eps_fd = Y_n[i_var] * Y_n[i_var];
        eps_fd = sqrt( eps_fd ) * 1e-4;

        // positive + eps
        for ( int j = 0; j < n_var; j++ )
            phi[j] = Y_n[j];

        phi[i_var] += 0.5 * eps_fd;
        phi[i_var]  = u_max( phi_bounds[i_var*BOUNDDIM], u_min( phi_bounds[i_var*BOUNDDIM+1], phi[i_var] ) );

        reactor_function_pointer( tpdt_loc );

        for ( int j = 0; j < n_var; j++ )
        {
            jac[i_var*n_var+j] = -phi_dt[j];
        }

        // negative + eps
        for ( int j = 0; j < n_var; j++ )
            phi[j] = Y_n[j];

        phi[i_var] -= 0.5 * eps_fd;
        phi[i_var]  = u_max( phi_bounds[i_var*BOUNDDIM], u_min( phi_bounds[i_var*BOUNDDIM+1], phi[i_var] ) );

        reactor_function_pointer( tpdt_loc );

        for ( int j = 0; j < n_var; j++ )
        {
            jac[i_var*n_var+j] += phi_dt[j];
            jac[i_var*n_var+j] *= bdf_b_loc / (eps_fd + SMALL);
        }

        jac[i_var*n_var+i_var] += 1.0 / dt_loc;
    }
}

int matrix_vector_numerical( double *x, double *b, int n_var, int m )
{
#ifdef DEBUG
    u_unused( m );
#endif /* DEBUG */

    for ( int j = 0; j < n_var; j++ )
    {
        b[j] = jac[j] * x[0];
    }

    for ( int i_var = 1; i_var < n_var; i_var++ )
    {
        for ( int j = 0; j < n_var; j++ )
        {
            b[j] += jac[i_var*n_var+j] * x[i_var];
        }
    }

    return 0;
}