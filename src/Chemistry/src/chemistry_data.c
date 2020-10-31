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
void deallocate_elements( Elements_t **elements );
void print_elements( Elements_t *elements );
void deallocate_specii( Specii_t **specii );
void print_specii( Specii_t *specii );
void deallocate_reactions( Reactions_t **reactions );
void print_reactions( Reactions_t *reactions );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
Chemistry_t *allocate_chemistry()
{
    Chemistry_t *tmp = allocate( sizeof( Chemistry_t ) );

    tmp->elements   = NULL;
    tmp->specii     = NULL;
    tmp->reactions  = NULL;

    return tmp;
}

void print_chemistry_info( Chemistry_t *chemistry )
{
    if (chemistry == NULL) return;
    printf_r( "chemistry\n" );

    print_elements( chemistry->elements );
    print_specii( chemistry->specii );
    print_reactions( chemistry->reactions );
}

void deallocate_chemistry( Chemistry_t **chemistry )
{
    if ((*chemistry) == NULL) return;

    deallocate_elements( &(*chemistry)->elements );
    deallocate_specii( &(*chemistry)->specii );
    deallocate_reactions( &(*chemistry)->reactions );

    deallocate( *chemistry );
}

Elements_t *allocate_elements( Chemistry_t *chemistry, int n_elements, int max_name_length )
{
    chemistry->elements = allocate( sizeof( Elements_t ) );
    Elements_t *elements = chemistry->elements;

    elements->n_elements    = n_elements;

    elements->symbol        = allocate_hdf5_string_buffer( n_elements, max_name_length );
    elements->mass          = allocate( sizeof( double ) * n_elements );

    return elements;
}

Specii_t *allocate_specii( Chemistry_t *chemistry, int n_specii, int max_name_length, int max_reac_points )
{
    chemistry->specii = allocate( sizeof( Specii_t ) );
    Specii_t *specii = chemistry->specii;

    Elements_t *elements    = chemistry->elements;
    int n_elements          = elements->n_elements;

    specii->n_specii        = n_specii;
    specii->max_reac_points = max_reac_points;

    specii->symbol          = allocate_hdf5_string_buffer( n_specii, max_name_length );
    specii->is_inert        = allocate( sizeof( int ) * n_specii );

    specii->composition     = allocate( sizeof( int ) * n_elements * n_specii );
    specii->phase           = allocate( sizeof( int ) * n_specii );
    specii->bounds          = allocate( sizeof( double ) * BOUNDS * n_specii );
    specii->coeff_high      = allocate( sizeof( double ) * NASA * n_specii );
    specii->coeff_low       = allocate( sizeof( double ) * NASA * n_specii );

    specii->geom            = allocate( sizeof( int ) * n_specii );
    specii->pot_lj          = allocate( sizeof( double ) * n_specii );
    specii->col_lj          = allocate( sizeof( double ) * n_specii );
    specii->dip_mo          = allocate( sizeof( double ) * n_specii );
    specii->pol             = allocate( sizeof( double ) * n_specii );
    specii->rot_rel         = allocate( sizeof( double ) * n_specii );

    specii->n_reac_points   = allocate( sizeof( int ) * n_specii );
    specii->reac_points     = allocate( sizeof( int ) * max_reac_points * n_specii );
    specii->nu_reac_points  = allocate( sizeof( double ) * max_reac_points * n_specii );

    set_value_int_n( 0, specii->n_reac_points, n_specii );

    specii->molar_mass      = allocate( sizeof( double ) * n_specii );
    specii->molecule_weight = allocate( sizeof( double ) * n_specii );
    specii->Rsp             = allocate( sizeof( double ) * n_specii );

    specii->omega           = allocate( sizeof( double ) * n_specii );

    return specii;
}

Reactions_t *allocate_reactions( Chemistry_t *chemistry, int n_reactions,
    int max_reactants, int max_products, int max_troe_coeff, int max_efficiencies )
{
    chemistry->reactions = allocate( sizeof( Reactions_t ) );
    Reactions_t *reactions = chemistry->reactions;

    reactions->n_reactions      = n_reactions;
    reactions->max_reactants    = max_reactants;
    reactions->max_products     = max_products;
    reactions->max_troe_coeff   = max_troe_coeff;
    reactions->max_efficiencies = max_efficiencies;

    reactions->type             = allocate( sizeof( int ) * n_reactions );
    reactions->falloff_species  = allocate( sizeof( int ) * n_reactions );

    reactions->arr_coeff        = allocate( sizeof( double ) * ARR * n_reactions );
    reactions->is_reversible    = allocate( sizeof( int ) * n_reactions );

    reactions->n_reactants      = allocate( sizeof( int ) * n_reactions );
    reactions->reactants        = allocate( sizeof( int ) * max_reactants * n_reactions );
    reactions->nu_reactants     = allocate( sizeof( double ) * max_reactants * n_reactions );
    reactions->ord_reactants    = allocate( sizeof( double ) * max_reactants * n_reactions );

    reactions->n_products       = allocate( sizeof( int ) * n_reactions );
    reactions->products         = allocate( sizeof( int ) * max_products * n_reactions );
    reactions->nu_products      = allocate( sizeof( double ) * max_products * n_reactions );
    reactions->ord_products     = allocate( sizeof( double ) * max_products * n_reactions );

    reactions->sum_nu           = allocate( sizeof( double ) * n_reactions );

    reactions->has_rev_arr      = allocate( sizeof( int ) * n_reactions );
    reactions->rev_arr_coeff    = allocate( sizeof( double ) * ARR * n_reactions );

    reactions->has_high_arr     = allocate( sizeof( int ) * n_reactions );
    reactions->has_low_arr      = allocate( sizeof( int ) * n_reactions );
    reactions->adv_arr_coeff    = allocate( sizeof( double ) * ARR * n_reactions );

    reactions->has_troe         = allocate( sizeof( int ) * n_reactions );
    reactions->n_troe_coeff     = allocate( sizeof( int ) * n_reactions );
    reactions->troe_coeff       = allocate( sizeof( double ) * max_troe_coeff * n_reactions );

    reactions->has_efficiencies = allocate( sizeof( int ) * n_reactions );
    reactions->n_efficiencies   = allocate( sizeof( int ) * n_reactions );
    reactions->sp_efficiencies  = allocate( sizeof( int ) * max_efficiencies * n_reactions );
    reactions->efficiencies     = allocate( sizeof( double ) * max_efficiencies * n_reactions );

    set_value_int_n( 0, reactions->n_reactants, n_reactions );
    set_value_int_n( 0, reactions->n_products, n_reactions );
    set_value_int_n( 0, reactions->n_efficiencies, n_reactions );

    reactions->n_reactions_three        = 0;
    reactions->idx_reactions_three      = NULL;
    reactions->n_reactions_low          = 0;
    reactions->idx_reactions_low        = NULL;
    reactions->n_reactions_high         = 0;
    reactions->idx_reactions_high       = NULL;
    reactions->n_reactions_troe         = 0;
    reactions->idx_reactions_troe       = NULL;
    reactions->n_reactions_rev          = 0;
    reactions->idx_reactions_rev        = NULL;
    reactions->n_reactions_rev_arr      = 0;
    reactions->idx_reactions_rev_arr    = NULL;

    reactions->q                        = allocate( sizeof( double ) * KINDIM * n_reactions );
    reactions->k                        = allocate( sizeof( double ) * KINDIM * n_reactions );
    reactions->pr                       = allocate( sizeof( double ) * n_reactions );

    return reactions;
}

void deallocate_elements( Elements_t **elements )
{
    if ((*elements) == NULL) return;

    deallocate_hdf5_string_buffer( &(*elements)->symbol );
    deallocate( (*elements)->mass );

    deallocate( *elements );
}

void print_elements( Elements_t *elements )
{
    if (elements == NULL) return;
    printf_r( "Elements\n" );

    printf_r( "n_elementss = %d\n", elements->n_elements );
}

void deallocate_specii( Specii_t **specii )
{
    if ((*specii) == NULL) return;

    deallocate_hdf5_string_buffer( &(*specii)->symbol );
    deallocate( (*specii)->is_inert );

    deallocate( (*specii)->composition );
    deallocate( (*specii)->phase );
    deallocate( (*specii)->bounds );
    deallocate( (*specii)->coeff_high );
    deallocate( (*specii)->coeff_low );

    deallocate( (*specii)->geom );
    deallocate( (*specii)->pot_lj );
    deallocate( (*specii)->col_lj );
    deallocate( (*specii)->dip_mo );
    deallocate( (*specii)->pol );
    deallocate( (*specii)->rot_rel );

    deallocate( (*specii)->n_reac_points );
    deallocate( (*specii)->reac_points );
    deallocate( (*specii)->nu_reac_points );

    deallocate( (*specii)->molar_mass );
    deallocate( (*specii)->molecule_weight );
    deallocate( (*specii)->Rsp );

    deallocate( (*specii)->omega );

    deallocate( *specii );
}

void print_specii( Specii_t *specii )
{
    if (specii == NULL) return;
    printf_r( "Specii\n" );

    printf_r( "n_specii = %d\n", specii->n_specii );
}

void deallocate_reactions( Reactions_t **reactions )
{
    if ((*reactions) == NULL) return;

    deallocate( (*reactions)->type );
    deallocate( (*reactions)->falloff_species );

    deallocate( (*reactions)->arr_coeff );
    deallocate( (*reactions)->is_reversible );

    deallocate( (*reactions)->n_reactants );
    deallocate( (*reactions)->reactants );
    deallocate( (*reactions)->nu_reactants );
    deallocate( (*reactions)->ord_reactants );

    deallocate( (*reactions)->n_products );
    deallocate( (*reactions)->products );
    deallocate( (*reactions)->nu_products );
    deallocate( (*reactions)->ord_products );

    deallocate( (*reactions)->sum_nu );

    deallocate( (*reactions)->has_rev_arr );
    deallocate( (*reactions)->rev_arr_coeff );

    deallocate( (*reactions)->has_high_arr );
    deallocate( (*reactions)->has_low_arr );
    deallocate( (*reactions)->adv_arr_coeff );

    deallocate( (*reactions)->has_troe );
    deallocate( (*reactions)->n_troe_coeff );
    deallocate( (*reactions)->troe_coeff );

    deallocate( (*reactions)->has_efficiencies );
    deallocate( (*reactions)->n_efficiencies );
    deallocate( (*reactions)->sp_efficiencies );
    deallocate( (*reactions)->efficiencies );

    deallocate( (*reactions)->idx_reactions_three );
    deallocate( (*reactions)->idx_reactions_low );
    deallocate( (*reactions)->idx_reactions_high );
    deallocate( (*reactions)->idx_reactions_troe );
    deallocate( (*reactions)->idx_reactions_rev );
    deallocate( (*reactions)->idx_reactions_rev_arr );

    deallocate( (*reactions)->q );
    deallocate( (*reactions)->k );
    deallocate( (*reactions)->pr );

    deallocate( *reactions );
}

void print_reactions( Reactions_t *reactions )
{
    if (reactions == NULL) return;
    printf_r( "Reactions\n" );

    printf_r( "n_reactions          = %d\n", reactions->n_reactions );
    printf_r( "n_reactions_three    = %d\n", reactions->n_reactions_three );
    printf_r( "n_reactions_low      = %d\n", reactions->n_reactions_low );
    printf_r( "n_reactions_high     = %d\n", reactions->n_reactions_high );
    printf_r( "n_reactions_troe     = %d\n", reactions->n_reactions_troe );
    printf_r( "n_reactions_rev      = %d\n", reactions->n_reactions_rev );
    printf_r( "n_reactions_rev_arr  = %d\n", reactions->n_reactions_rev_arr );
}