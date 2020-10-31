//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <math.h>
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
void calc_rate_of_progress( double *C, double T, Chemistry_t *chemistry );
void calc_reaction_rates( double T, Chemistry_t *chemistry );
double calc_arr( double T, double *coeff );
double calc_troe( double T, double pr, double *coeff );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void calc_production_rate( double *C, double T, Chemistry_t *chemistry )
{
    calc_rate_of_progress( C, T, chemistry );

    Specii_t *specii        = chemistry->specii;
    Reactions_t *reactions  = chemistry->reactions;
    int n_specii            = specii->n_specii;
    int max_reac_points     = specii->max_reac_points;

    for ( int i = 0; i < n_specii; i++ )
    {
        specii->omega[i] = 0.0;

        for ( int j = 0; j < specii->n_reac_points[i]; j++ )
        {
            int idx_r = specii->reac_points[i*max_reac_points+j] * KINDIM;

            specii->omega[i] += specii->nu_reac_points[i*max_reac_points+j] *
                (reactions->q[idx_r] - reactions->q[idx_r+1]);
        }

        specii->omega[i] *= specii->molar_mass[i];
    }
}

void calc_rate_of_progress( double *C, double T, Chemistry_t *chemistry )
{
    calc_reaction_rates( T, chemistry );

    Specii_t *specii        = chemistry->specii;
    Reactions_t *reactions  = chemistry->reactions;
    int n_specii            = specii->n_specii;
    int n_reactions         = reactions->n_reactions;
    int n_reactions_rev     = reactions->n_reactions_rev;
    int n_reactions_rev_arr = reactions->n_reactions_rev_arr;
    int n_reactions_three   = reactions->n_reactions_three;
    int n_reactions_low     = reactions->n_reactions_low;
    int n_reactions_high    = reactions->n_reactions_high;
    int n_reactions_troe    = reactions->n_reactions_troe;
    int max_reactants       = reactions->max_reactants;
    int max_products        = reactions->max_products;
    int max_troe_coeff      = reactions->max_troe_coeff;
    int max_efficiencies    = reactions->max_efficiencies;

    double C_sum = sum_n( C, n_specii );
    int use_sp[n_specii];

    for ( int i = 0; i < n_specii; i++ )
        use_sp[i] = (C[i] > YSMALL);

    // assemble the forward reaction rate or progress
    for ( int i = 0; i < n_reactions; i++ )
    {
        int is_zero = 0;

        double tmp = 0.0;
        for ( int j = 0; j < reactions->n_reactants[i]; j++ )
        {
            int idx_sp = i*max_reactants+j;
            if (use_sp[reactions->reactants[idx_sp]] == 0) is_zero = 1;
            tmp += log( C[reactions->reactants[idx_sp]] ) * reactions->ord_reactants[idx_sp];
        }

        reactions->q[i*KINDIM] = (is_zero == 0) ? reactions->k[i*KINDIM] * exp( tmp ) : 0.0;
    }

    // assemble the backward reaction rate or progress
    for ( int ii = 0; ii < n_reactions_rev; ii++ )
    {
        int i       = reactions->idx_reactions_rev[i];
        int is_zero = 0;
        double tmp  = 0.0;

        for ( int j = 0; j < reactions->n_products[i]; j++ )
        {
            int idx_sp = i*max_products+j;
            if (use_sp[reactions->products[idx_sp]] == 0) is_zero = 1;
            tmp += log( C[reactions->products[idx_sp]] ) * reactions->ord_products[idx_sp];
        }

        reactions->q[i*KINDIM+1] = (is_zero == 0) ? reactions->k[i*KINDIM+1] * exp( tmp ) : 0.0;
    }

    for ( int ii = 0; ii < n_reactions_rev_arr; ii++ )
    {
        int i       = reactions->idx_reactions_rev_arr[i];
        int is_zero = 0;
        double tmp  = 0.0;

        for ( int j = 0; j < reactions->n_products[i]; j++ )
        {
            int idx_sp = i*max_products+j;
            if (use_sp[reactions->products[idx_sp]] == 0) is_zero = 1;
            tmp += log( C[reactions->products[idx_sp]] ) * reactions->ord_products[idx_sp];
        }

        reactions->q[i*KINDIM+1] = (is_zero == 0) ? reactions->k[i*KINDIM+1] * exp( tmp ) : 0.0;
    }

    // three-body type reactions
    for ( int ii = 0; ii < n_reactions_three; ii++ )
    {
        int i = reactions->idx_reactions_three[i];

        C[n_specii] = C_sum;
        for ( int j = 0; j < reactions->n_efficiencies[i]; j++ )
        {
            int idx_sp = i*max_efficiencies+j;
            C[n_specii] += (reactions->efficiencies[idx_sp] - 1.0) * C[reactions->sp_efficiencies[idx_sp]];
        }

        reactions->q[i*KINDIM]      *= C[reactions->falloff_species[i]];
        reactions->q[i*KINDIM+1]    *= C[reactions->falloff_species[i]];
    }

    // pressure type reactions (unimolecular/recombination fall-off reactions)
    for ( int ii = 0; ii < n_reactions_low; ii++ )
    {
        int i = reactions->idx_reactions_low[i];

        C[n_specii] = C_sum;
        for ( int j = 0; j < reactions->n_efficiencies[i]; j++ )
        {
            int idx_sp = i*max_efficiencies+j;
            C[n_specii] += (reactions->efficiencies[idx_sp] - 1.0) * C[reactions->sp_efficiencies[idx_sp]];
        }

        reactions->pr[i] = C[reactions->falloff_species[i]] *
            calc_arr( T, &reactions->adv_arr_coeff[i*ARR] ) / reactions->k[i*KINDIM];

        double tmp                   = reactions->pr[i] / (1.0 + reactions->pr[i]);
        reactions->q[i*KINDIM]      *= tmp;
        reactions->q[i*KINDIM+1]    *= tmp;
    }

    // pressure type reactions (chemically activated bimolecular reactions)
    for ( int ii = 0; ii < n_reactions_high; ii++ )
    {
        int i = reactions->idx_reactions_high[i];

        C[n_specii] = C_sum;
        for ( int j = 0; j < reactions->n_efficiencies[i]; j++ )
        {
            int idx_sp = i*max_efficiencies+j;
            C[n_specii] += (reactions->efficiencies[idx_sp] - 1.0) * C[reactions->sp_efficiencies[idx_sp]];
        }

        reactions->pr[i] = C[reactions->falloff_species[i]] * reactions->k[i*KINDIM] /
            calc_arr( T, &reactions->adv_arr_coeff[i*ARR] );

        double tmp                   = 1.0 / (1.0 + reactions->pr[i]);
        reactions->q[i*KINDIM]      *= tmp;
        reactions->q[i*KINDIM+1]    *= tmp;
    }

    // pressure type reactions (Troe)
    for ( int ii = 0; ii < n_reactions_troe; ii++ )
    {
        int i       = reactions->idx_reactions_troe[i];
        double tmp  = calc_troe( T, reactions->pr[i], &reactions->troe_coeff[i*max_troe_coeff] );

        reactions->q[i*KINDIM]      *= tmp;
        reactions->q[i*KINDIM+1]    *= tmp;
    }
}

void calc_reaction_rates( double T, Chemistry_t *chemistry )
{
    Specii_t *specii        = chemistry->specii;
    Reactions_t *reactions  = chemistry->reactions;
    int n_specii            = specii->n_specii;
    int n_reactions         = reactions->n_reactions;
    int n_reactions_rev     = reactions->n_reactions_rev;
    int n_reactions_rev_arr = reactions->n_reactions_rev_arr;
    int max_reactants       = reactions->max_reactants;
    int max_products        = reactions->max_products;

    double s_prmt = P0 / (RM * T);
    double g_rt[n_specii];

    for ( int i = 0; i < n_specii; i++ )
    {
        g_rt[i] = calc_sp_g_rt( i, T, chemistry );
    }

    // forward reaction rate
    for ( int i = 0; i < n_reactions; i++ )
    {
        reactions->k[i*KINDIM] = calc_arr( T, &reactions->arr_coeff[i*ARR] );
    }

    // backward reaction rate from equilibrium constant
    for ( int ii = 0; ii < n_reactions_rev; ii++ )
    {
        int i = reactions->idx_reactions_rev[ii];

        double tmp = 0.0;
        for ( int j = 0; j < reactions->n_reactants[i]; j++ )
        {
            int idx_sp = i*max_reactants+j;
            tmp -= reactions->nu_reactants[idx_sp] * g_rt[idx_sp];
        }

        for ( int j = 0; j < reactions->n_products[i]; j++ )
        {
            int idx_sp = i*max_products+j;
            tmp += reactions->nu_products[idx_sp] * g_rt[idx_sp];
        }

        tmp = exp( tmp );
        reactions->k[i*KINDIM+1] = reactions->k[i*KINDIM] / (tmp * pow( s_prmt, reactions->sum_nu[i] ));
    }

    // backward reaction rate from reverse rate constant
    for ( int ii = 0; ii < n_reactions_rev_arr; ii++ )
    {
        int i = reactions->idx_reactions_rev_arr[ii];
        reactions->k[i*KINDIM+1] = calc_arr( T, &reactions->rev_arr_coeff[i*ARR] );
    }
}

double calc_arr( double T, double *coeff )
{
    return coeff[0] * exp( coeff[1] * log( T ) - coeff[2] / (RM * T) );
}

double calc_troe( double T, double pr, double *coeff )
{
    double l10_pr       = log10( pr );
    double fcn_tr       = (1. - coeff[0]) * exp( -T / coeff[1] ) + coeff[0] * exp( -T / coeff[2] ) + exp( -coeff[3] / T );
    double l10_fcn_tr   = log10( fcn_tr );

    double c_tr         = -0.4 - 0.67 * l10_fcn_tr;
    double nk_tr        = 0.75 - 1.27 * l10_fcn_tr;
    double d_tr         = 0.14;

    double bk_tr        = l10_pr + c_tr;
    double a_tr         = bk_tr /(nk_tr - d_tr * (l10_pr + c_tr));

    return pow( 10, l10_fcn_tr / (1. + a_tr * a_tr) );
}