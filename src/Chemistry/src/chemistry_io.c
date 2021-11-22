/*******************************************************************************
 * @file chemistry.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "chemistry_private.h"

/*******************************************************************************
 * @brief Read chemistry file and fill chemistry structure
 * @param chemistry_file
 * @return chemistry_t*
 ******************************************************************************/
chemistry_t *read_chemistry_data(cstring_t chemistry_file)
{
    chemistry_t *chemistry = allocate_chemistry();

    hid_t file_id = open_hdf5_file(chemistry_file);

    /* the elements */
    {
        hsize_t dims[2];

        int n_elements = 0;

        hid_t group_id = open_hdf5_group(file_id, "ELEMENTS");
        GET_HDF5_ATTRIBUTE(group_id, "n_elements", HDF5Int, &n_elements);

        int max_length = get_hdf5_dataset_size(group_id, "symbol");
        elements_t *elements =
            allocate_elements(chemistry, n_elements, max_length);

        dims[0] = elements->n_elements;
        GET_HDF5_DATASET_N(group_id, "symbol", HDF5String,
                           dims[0], elements->symbol);

        dims[0] = elements->n_elements;
        GET_HDF5_DATASET_N(group_id, "mass", HDF5Double,
                           dims[0], elements->mass);

        close_hdf5_group(group_id);
    }

    /* the specii */
    {
        hsize_t dims[2];

        int n_specii = 0;
        int max_reac_points = 0;

        hid_t group_id = open_hdf5_group(file_id, "SPECII");
        GET_HDF5_ATTRIBUTE(group_id, "n_specii", HDF5Int, &n_specii);
        GET_HDF5_ATTRIBUTE(group_id, "max_reac_points", HDF5Int,
                           &max_reac_points);

        int max_length = get_hdf5_dataset_size(group_id, "symbol");
        specii_t *specii =
            allocate_specii(chemistry, n_specii, max_length, max_reac_points);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "symbol", HDF5String,
                           dims[0], specii->symbol);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "is_inert", HDF5Int,
                           dims[0], specii->is_inert);

        dims[0] = specii->n_specii;
        dims[1] = chemistry->elements->n_elements;
        GET_HDF5_DATASET_N_M(group_id, "composition", HDF5Int,
                             dims, specii->composition);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "phase", HDF5Int,
                           dims[0], specii->phase);

        dims[0] = specii->n_specii;
        dims[1] = BOUNDS;
        GET_HDF5_DATASET_N_M(group_id, "bounds", HDF5Double,
                             dims, specii->bounds);

        dims[0] = specii->n_specii;
        dims[1] = NASA;
        GET_HDF5_DATASET_N_M(group_id, "coeff_high", HDF5Double,
                             dims, specii->coeff_high);

        dims[0] = specii->n_specii;
        dims[1] = NASA;
        GET_HDF5_DATASET_N_M(group_id, "coeff_low", HDF5Double,
                             dims, specii->coeff_low);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "geom", HDF5Int, dims[0], specii->geom);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "pot_lj", HDF5Double,
                           dims[0], specii->pot_lj);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "col_lj", HDF5Double,
                           dims[0], specii->col_lj);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "dip_mo", HDF5Double,
                           dims[0], specii->dip_mo);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "pol", HDF5Double, dims[0], specii->pol);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "rot_rel", HDF5Double,
                           dims[0], specii->rot_rel);

        dims[0] = specii->n_specii;
        GET_HDF5_DATASET_N(group_id, "n_reac_points", HDF5Int,
                           dims[0], specii->n_reac_points);

        dims[0] = specii->n_specii;
        dims[1] = specii->max_reac_points;
        GET_HDF5_DATASET_N_M(group_id, "reac_points", HDF5Int,
                             dims, specii->reac_points);

        dims[0] = specii->n_specii;
        dims[1] = specii->max_reac_points;
        GET_HDF5_DATASET_N_M(group_id, "nu_reac_points", HDF5Double,
                             dims, specii->nu_reac_points);

        close_hdf5_group(group_id);
    }

    /* the reactions */
    {
        hsize_t dims[2];

        int n_reactions = 0;
        int max_reactants = 0;
        int max_products = 0;
        int max_troe_coeff = 0;
        int max_efficiencies = 0;

        hid_t group_id = open_hdf5_group(file_id, "REACTIONS");
        GET_HDF5_ATTRIBUTE(group_id, "n_reactions", HDF5Int, &n_reactions);
        GET_HDF5_ATTRIBUTE(group_id, "max_reactants", HDF5Int, &max_reactants);
        GET_HDF5_ATTRIBUTE(group_id, "max_products", HDF5Int, &max_products);
        GET_HDF5_ATTRIBUTE(group_id, "max_troe_coeff", HDF5Int,
                           &max_troe_coeff);
        GET_HDF5_ATTRIBUTE(group_id, "max_efficiencies", HDF5Int,
                           &max_efficiencies);

        reactions_t *reactions =
            allocate_reactions(chemistry, n_reactions,
                               max_reactants, max_products,
                               max_troe_coeff, max_efficiencies);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "type", HDF5Int, dims[0], reactions->type);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "falloff_species", HDF5Int,
                           dims[0], reactions->falloff_species);

        dims[0] = reactions->n_reactions;
        dims[1] = ARR;
        GET_HDF5_DATASET_N_M(group_id, "arr_coeff", HDF5Double,
                             dims, reactions->arr_coeff);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "is_reversible", HDF5Int,
                           dims[0], reactions->is_reversible);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "n_reactants", HDF5Int,
                           dims[0], reactions->n_reactants);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_reactants;
        GET_HDF5_DATASET_N_M(group_id, "reactants", HDF5Int,
                             dims, reactions->reactants);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_reactants;
        GET_HDF5_DATASET_N_M(group_id, "nu_reactants", HDF5Double,
                             dims, reactions->nu_reactants);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_reactants;
        GET_HDF5_DATASET_N_M(group_id, "ord_reactants", HDF5Double,
                             dims, reactions->ord_reactants);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "n_products", HDF5Int,
                           dims[0], reactions->n_products);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_products;
        GET_HDF5_DATASET_N_M(group_id, "products", HDF5Int,
                             dims, reactions->products);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_reactants;
        GET_HDF5_DATASET_N_M(group_id, "nu_products", HDF5Double,
                             dims, reactions->nu_products);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_reactants;
        GET_HDF5_DATASET_N_M(group_id, "ord_products", HDF5Double,
                             dims, reactions->ord_products);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "sum_nu", HDF5Double,
                           dims[0], reactions->sum_nu);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "has_rev_arr", HDF5Int,
                           dims[0], reactions->has_rev_arr);

        dims[0] = reactions->n_reactions;
        dims[1] = ARR;
        GET_HDF5_DATASET_N_M(group_id, "rev_arr_coeff", HDF5Double,
                             dims, reactions->rev_arr_coeff);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "has_high_arr", HDF5Int,
                           dims[0], reactions->has_high_arr);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "has_low_arr", HDF5Int,
                           dims[0], reactions->has_low_arr);

        dims[0] = reactions->n_reactions;
        dims[1] = ARR;
        GET_HDF5_DATASET_N_M(group_id, "adv_arr_coeff", HDF5Double,
                             dims, reactions->adv_arr_coeff);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "has_troe", HDF5Int,
                           dims[0], reactions->has_troe);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "n_troe_coeff", HDF5Int,
                           dims[0], reactions->n_troe_coeff);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_troe_coeff;
        GET_HDF5_DATASET_N_M(group_id, "troe_coeff", HDF5Double,
                             dims, reactions->troe_coeff);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "has_efficiencies", HDF5Int,
                           dims[0], reactions->has_efficiencies);

        dims[0] = reactions->n_reactions;
        GET_HDF5_DATASET_N(group_id, "n_efficiencies", HDF5Int,
                           dims[0], reactions->n_efficiencies);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_efficiencies;
        GET_HDF5_DATASET_N_M(group_id, "sp_efficiencies", HDF5Int,
                             dims, reactions->sp_efficiencies);

        dims[0] = reactions->n_reactions;
        dims[1] = reactions->max_efficiencies;
        GET_HDF5_DATASET_N_M(group_id, "efficiencies", HDF5Double,
                             dims, reactions->efficiencies);

        close_hdf5_group(group_id);
    }

    close_hdf5_file(file_id);

    complete_chemistry(chemistry);
    return chemistry;
}
