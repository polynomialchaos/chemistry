//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "chemistry_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
double calc_ig_p( double rho, double R, double T );
double calc_ig_rho( double p, double R, double T );
double calc_ig_T( double p, double rho, double R );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
State_t *allocate_state( Chemistry_t *chemistry )
{
    State_t *state      = allocate( sizeof( State_t ) );
    state->chemistry    = chemistry;

    state->Y            = allocate( sizeof( double ) * chemistry->specii->n_specii );
    state->C            = allocate( sizeof( double ) * (chemistry->specii->n_specii + 1) ); // 3rd body concentration
    state->X            = allocate( sizeof( double ) * chemistry->specii->n_specii );

    set_value_n( 0.0, state->Y, chemistry->specii->n_specii );
    set_value_n( 0.0, state->C, (chemistry->specii->n_specii + 1) );
    set_value_n( 0.0, state->X, chemistry->specii->n_specii );

    return state;
}

void deallocate_state( State_t **state )
{
    if ((*state) == NULL) return;

    deallocate( (*state)->Y );
    deallocate( (*state)->C );
    deallocate( (*state)->X );

    deallocate( (*state) );
}

void update_state_isochoric( double p, double rho, double T, double *Y, State_t *state )
{
#ifdef DEBUG
    u_unused( p );
#endif /* DEBUG */

    Chemistry_t *chemistry  = state->chemistry;

    state->rho  = rho;
    state->T    = T;

    correct_fraction( Y, state->Y, chemistry );

    state->R    = calc_mix_R( state->Y, chemistry );
    state->p    = calc_ig_p( state->rho, state->R, state->T );

    conv_mass_to_conc( state->Y, state->rho, state->C, chemistry );

    state->molar_mass   = calc_mix_molar_mass_from_Y( state->Y, chemistry );

    conv_mass_to_mole( state->Y, state->molar_mass, state->X, chemistry );

    state->cp   = calc_mix_cp( state->Y, state->T, chemistry );
    state->cv   = calc_mix_cv( state->cp, state->R );
    state->h    = calc_mix_h( state->Y, state->T, chemistry );
    state->s    = calc_mix_s( state->X, state->Y, state->p, state->T, chemistry );
    state->u    = calc_mix_u( state->h, state->R, state->T );
    state->g    = calc_mix_g( state->h, state->s, state->T );
}

void update_state_isobaric( double p, double rho, double T, double *Y, State_t *state )
{
#ifdef DEBUG
    u_unused( rho );
#endif /* DEBUG */

    Chemistry_t *chemistry  = state->chemistry;

    state->p    = p;
    state->T    = T;

    correct_fraction( Y, state->Y, chemistry );

    state->R    = calc_mix_R( state->Y, chemistry );
    state->rho  = calc_ig_rho( state->p, state->R, state->T );

    conv_mass_to_conc( state->Y, state->rho, state->C, chemistry );

    state->molar_mass   = calc_mix_molar_mass_from_Y( state->Y, chemistry );

    conv_mass_to_mole( state->Y, state->molar_mass, state->X, chemistry );

    state->cp   = calc_mix_cp( state->Y, state->T, chemistry );
    state->cv   = calc_mix_cv( state->cp, state->R );
    state->h    = calc_mix_h( state->Y, state->T, chemistry );
    state->s    = calc_mix_s( state->X, state->Y, state->p, state->T, chemistry );
    state->u    = calc_mix_u( state->h, state->R, state->T );
    state->g    = calc_mix_g( state->h, state->s, state->T );
}

double calc_ig_p( double rho, double R, double T )
{
    return rho * R * T;
}

double calc_ig_rho( double p, double R, double T )
{
    return p / (R * T);
}

double calc_ig_T( double p, double rho, double R )
{
    return p / (rho * R);
}