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
enum ReactionType {
    ReactionDefault, ReactionThreeBody, ReactionPressure, ReactionTypeMax
};

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void calc_chemistry_data( Chemistry_t *chemistry );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
Chemistry_t *read_chemistry_data( const_string_t chemistry_file )
{
    Chemistry_t *chemistry = allocate_chemistry();

    hid_t file_id = open_hdf5_file( chemistry_file );

        // the elements
        {
            hsize_t dims[2];

            int n_elements = 0;

            hid_t group_id = open_hdf5_group( file_id, "ELEMENTS" );
                get_hdf5_attribute( group_id, "n_elements", HDF5Int, &n_elements );

                int max_length = get_hdf5_dataset_size( group_id, "symbol" );
                Elements_t *elements = allocate_elements( chemistry, n_elements, max_length );

                dims[0] = elements->n_elements;
                get_hdf5_dataset_n( group_id, "symbol", HDF5String, elements->symbol[0], 1, dims );

                dims[0] = elements->n_elements;
                get_hdf5_dataset_n( group_id, "mass", HDF5Double, elements->mass, 1, dims );
            close_hdf5_group( group_id );
        }

        // the specii
        {
            hsize_t dims[2];

            int n_specii = 0;
            int max_reac_points = 0;

            hid_t group_id = open_hdf5_group( file_id, "SPECII" );
                get_hdf5_attribute( group_id, "n_specii", HDF5Int, &n_specii );
                get_hdf5_attribute( group_id, "max_reac_points", HDF5Int, &max_reac_points );

                int max_length = get_hdf5_dataset_size( group_id, "symbol" );
                Specii_t *specii = allocate_specii( chemistry, n_specii, max_length, max_reac_points );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "symbol", HDF5String, specii->symbol[0], 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "is_inert", HDF5Int, specii->is_inert, 1, dims );

                dims[0] = specii->n_specii;
                dims[1] = chemistry->elements->n_elements;
                get_hdf5_dataset_n_m( group_id, "composition", HDF5Int, specii->composition, 2, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "phase", HDF5Int, specii->phase, 1, dims );

                dims[0] = specii->n_specii;
                dims[1] = BOUNDS;
                get_hdf5_dataset_n_m( group_id, "bounds", HDF5Double, specii->bounds, 2, dims );

                dims[0] = specii->n_specii;
                dims[1] = NASA;
                get_hdf5_dataset_n_m( group_id, "coeff_high", HDF5Double, specii->coeff_high, 2, dims );

                dims[0] = specii->n_specii;
                dims[1] = NASA;
                get_hdf5_dataset_n_m( group_id, "coeff_low", HDF5Double, specii->coeff_low, 2, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "geom", HDF5Int, specii->geom, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "pot_lj", HDF5Double, specii->pot_lj, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "col_lj", HDF5Double, specii->col_lj, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "dip_mo", HDF5Double, specii->dip_mo, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "pol", HDF5Double, specii->pol, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "rot_rel", HDF5Double, specii->rot_rel, 1, dims );

                dims[0] = specii->n_specii;
                get_hdf5_dataset_n( group_id, "n_reac_points", HDF5Int, specii->n_reac_points, 1, dims );

                dims[0] = specii->n_specii;
                dims[1] = specii->max_reac_points;
                get_hdf5_dataset_n_m( group_id, "reac_points", HDF5Int, specii->reac_points, 2, dims );

                dims[0] = specii->n_specii;
                dims[1] = specii->max_reac_points;
                get_hdf5_dataset_n_m( group_id, "nu_reac_points", HDF5Double, specii->nu_reac_points, 2, dims );
            close_hdf5_group( group_id );
        }

        // the reactions
        {
            hsize_t dims[2];

            int n_reactions = 0;
            int max_reactants = 0;
            int max_products = 0;
            int max_troe_coeff = 0;
            int max_efficiencies = 0;

            hid_t group_id = open_hdf5_group( file_id, "REACTIONS" );
                get_hdf5_attribute( group_id, "n_reactions", HDF5Int, &n_reactions );
                get_hdf5_attribute( group_id, "max_reactants", HDF5Int, &max_reactants );
                get_hdf5_attribute( group_id, "max_products", HDF5Int, &max_products );
                get_hdf5_attribute( group_id, "max_troe_coeff", HDF5Int, &max_troe_coeff );
                get_hdf5_attribute( group_id, "max_efficiencies", HDF5Int, &max_efficiencies );

                Reactions_t *reactions = allocate_reactions( chemistry, n_reactions,
                    max_reactants, max_products, max_troe_coeff, max_efficiencies );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "type", HDF5Int, reactions->type, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "falloff_species", HDF5Int, reactions->falloff_species, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = ARR;
                get_hdf5_dataset_n_m( group_id, "arr_coeff", HDF5Double, reactions->arr_coeff, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "is_reversible", HDF5Int, reactions->is_reversible, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "n_reactants", HDF5Int, reactions->n_reactants, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_reactants;
                get_hdf5_dataset_n_m( group_id, "reactants", HDF5Int, reactions->reactants, 2, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_reactants;
                get_hdf5_dataset_n_m( group_id, "nu_reactants", HDF5Double, reactions->nu_reactants, 2, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_reactants;
                get_hdf5_dataset_n_m( group_id, "ord_reactants", HDF5Double, reactions->ord_reactants, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "n_products", HDF5Int, reactions->n_products, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_products;
                get_hdf5_dataset_n_m( group_id, "products", HDF5Int, reactions->products, 2, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_reactants;
                get_hdf5_dataset_n_m( group_id, "nu_products", HDF5Double, reactions->nu_products, 2, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_reactants;
                get_hdf5_dataset_n_m( group_id, "ord_products", HDF5Double, reactions->ord_products, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "sum_nu", HDF5Double, reactions->sum_nu, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "has_rev_arr", HDF5Int, reactions->has_rev_arr, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = ARR;
                get_hdf5_dataset_n_m( group_id, "rev_arr_coeff", HDF5Double, reactions->rev_arr_coeff, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "has_high_arr", HDF5Int, reactions->has_high_arr, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "has_low_arr", HDF5Int, reactions->has_low_arr, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = ARR;
                get_hdf5_dataset_n_m( group_id, "adv_arr_coeff", HDF5Double, reactions->adv_arr_coeff, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "has_troe", HDF5Int, reactions->has_troe, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "n_troe_coeff", HDF5Int, reactions->n_troe_coeff, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_troe_coeff;
                get_hdf5_dataset_n_m( group_id, "troe_coeff", HDF5Double, reactions->troe_coeff, 2, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "has_efficiencies", HDF5Int, reactions->has_efficiencies, 1, dims );

                dims[0] = reactions->n_reactions;
                get_hdf5_dataset_n( group_id, "n_efficiencies", HDF5Int, reactions->n_efficiencies, 1, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_efficiencies;
                get_hdf5_dataset_n_m( group_id, "sp_efficiencies", HDF5Int, reactions->sp_efficiencies, 2, dims );

                dims[0] = reactions->n_reactions;
                dims[1] = reactions->max_efficiencies;
                get_hdf5_dataset_n_m( group_id, "efficiencies", HDF5Double, reactions->efficiencies, 2, dims );
            close_hdf5_group( group_id );
        }

    close_hdf5_file( file_id );

    calc_chemistry_data( chemistry );
    return chemistry;
}

void calc_chemistry_data( Chemistry_t *chemistry )
{
    Elements_t *elements    = chemistry->elements;
    Specii_t *specii        = chemistry->specii;
    Reactions_t *reactions  = chemistry->reactions;
    int n_elements          = elements->n_elements;
    int n_specii            = specii->n_specii;
    int n_reactions         = reactions->n_reactions;

    for ( int i = 0; i < n_specii; i++ )
    {
        specii->molar_mass[i] = 0.0;
        for ( int j = 0; j < n_elements; j++ )
            specii->molar_mass[i] += specii->composition[i*n_elements+j] * elements->mass[j];

        specii->molecule_weight[i]  = specii->molar_mass[i] / NA;
        specii->Rsp[i]              = RM / specii->molar_mass[i];
    }

    if (n_reactions > 0)
    {
        reactions->n_reactions_three = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if (reactions->type[i] == ReactionThreeBody)
            {
                reactions->n_reactions_three += 1;

                reactions->idx_reactions_three = reallocate(
                    reactions->idx_reactions_three, sizeof( int ) * reactions->n_reactions_three );

                reactions->idx_reactions_three[reactions->n_reactions_three-1] = i;
            }

        reactions->n_reactions_low = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if ((reactions->type[i] == ReactionPressure) && (reactions->has_low_arr[i] == 1))
            {
                reactions->n_reactions_low += 1;

                reactions->idx_reactions_low = reallocate(
                    reactions->idx_reactions_low, sizeof( int ) * reactions->n_reactions_low );

                reactions->idx_reactions_low[reactions->n_reactions_low-1] = i;
            }

        reactions->n_reactions_high = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if ((reactions->type[i] == ReactionPressure) && (reactions->has_high_arr[i] == 1))
            {
                reactions->n_reactions_high += 1;

                reactions->idx_reactions_high = reallocate(
                    reactions->idx_reactions_high, sizeof( int ) * reactions->n_reactions_high );

                reactions->idx_reactions_high[reactions->n_reactions_high-1] = i;
            }

        reactions->n_reactions_troe = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if ((reactions->type[i] == ReactionPressure) && (reactions->has_troe[i] == 1))
            {
                reactions->n_reactions_troe += 1;

                reactions->idx_reactions_troe = reallocate(
                    reactions->idx_reactions_troe, sizeof( int ) * reactions->n_reactions_troe );

                reactions->idx_reactions_troe[reactions->n_reactions_troe-1] = i;
            }

        reactions->n_reactions_rev = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if ((reactions->is_reversible[i] == 1) && (reactions->has_rev_arr[i] == 0))
            {
                reactions->n_reactions_rev += 1;

                reactions->idx_reactions_rev = reallocate(
                    reactions->idx_reactions_rev, sizeof( int ) * reactions->n_reactions_rev );

                reactions->idx_reactions_rev[reactions->n_reactions_rev-1] = i;
            }

        reactions->n_reactions_rev_arr = 0;
        for ( int i = 0; i < n_reactions; i++ )
            if ((reactions->is_reversible[i] == 1) && (reactions->has_rev_arr[i] == 1))
            {
                reactions->n_reactions_rev_arr += 1;

                reactions->idx_reactions_rev_arr = reallocate(
                    reactions->idx_reactions_rev_arr, sizeof( int ) * reactions->n_reactions_rev_arr );

                reactions->idx_reactions_rev_arr[reactions->n_reactions_rev_arr-1] = i;
            }
    }
}

int get_species_index( Chemistry_t *chemistry, const_string_t symbol )
{
    Specii_t *specii    = chemistry->specii;
    int n_specii        = specii->n_specii;
    string_t *symbols   = specii->symbol;

    for ( int i = 0; i < n_specii; i++ )
    {
        if (is_equal( symbol, symbols[i] )) return i;
    }

    check_error( 0 );
    return -1;
}