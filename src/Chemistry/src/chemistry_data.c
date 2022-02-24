/*******************************************************************************
 * @file chemistry_data.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "chemistry_private.h"

#define PCSS '-'   /** Chemistry print separator symbol */
#define PCIN " * " /** Chemistry print intent */

/*******************************************************************************
 * @brief Allocate chemistry structure
 * @return chemistry_t*
 ******************************************************************************/
chemistry_t *allocate_chemistry()
{
    chemistry_t *tmp = BM_ALLOCATE(sizeof(chemistry_t));

    tmp->elements = NULL;
    tmp->specii = NULL;
    tmp->reactions = NULL;

    return tmp;
}

/*******************************************************************************
 * @brief Allocate elements structure
 * @param chemistry
 * @param n_elements
 * @param max_name_length
 * @return elements_t*
 ******************************************************************************/
elements_t *allocate_elements(chemistry_t *chemistry,
                              int n_elements, int max_name_length)
{
    chemistry->elements = BM_ALLOCATE(sizeof(elements_t));
    elements_t *elements = chemistry->elements;

    elements->n_elements = n_elements;

    elements->symbol =
        allocate_hdf5_string_buffer(n_elements, max_name_length, NULL);
    elements->mass = BM_ALLOCATE(sizeof(double) * n_elements);

    return elements;
}

/*******************************************************************************
 * @brief Allocate reactions structure
 * @param chemistry
 * @param n_reactions
 * @param max_reactants
 * @param max_products
 * @param max_troe_coeff
 * @param max_efficiencies
 * @return reactions_t*
 ******************************************************************************/
reactions_t *allocate_reactions(chemistry_t *chemistry, int n_reactions,
                                int max_reactants, int max_products,
                                int max_troe_coeff, int max_efficiencies)
{
    chemistry->reactions = BM_ALLOCATE(sizeof(reactions_t));
    reactions_t *reactions = chemistry->reactions;

    reactions->n_reactions = n_reactions;
    reactions->max_reactants = max_reactants;
    reactions->max_products = max_products;
    reactions->max_troe_coeff = max_troe_coeff;
    reactions->max_efficiencies = max_efficiencies;

    reactions->type = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->falloff_species = BM_ALLOCATE(sizeof(int) * n_reactions);

    reactions->arr_coeff = BM_ALLOCATE(sizeof(double) * ARR * n_reactions);
    reactions->is_reversible = BM_ALLOCATE(sizeof(int) * n_reactions);

    reactions->n_reactants = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->reactants = BM_ALLOCATE(sizeof(int) * max_reactants * n_reactions);
    reactions->nu_reactants =
        BM_ALLOCATE(sizeof(double) * max_reactants * n_reactions);
    reactions->ord_reactants =
        BM_ALLOCATE(sizeof(double) * max_reactants * n_reactions);

    reactions->n_products = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->products = BM_ALLOCATE(sizeof(int) * max_products * n_reactions);
    reactions->nu_products =
        BM_ALLOCATE(sizeof(double) * max_products * n_reactions);
    reactions->ord_products =
        BM_ALLOCATE(sizeof(double) * max_products * n_reactions);

    reactions->sum_nu = BM_ALLOCATE(sizeof(double) * n_reactions);

    reactions->has_rev_arr = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->rev_arr_coeff = BM_ALLOCATE(sizeof(double) * ARR * n_reactions);

    reactions->has_high_arr = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->has_low_arr = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->adv_arr_coeff = BM_ALLOCATE(sizeof(double) * ARR * n_reactions);

    reactions->has_troe = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->n_troe_coeff = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->troe_coeff =
        BM_ALLOCATE(sizeof(double) * max_troe_coeff * n_reactions);

    reactions->has_efficiencies = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->n_efficiencies = BM_ALLOCATE(sizeof(int) * n_reactions);
    reactions->sp_efficiencies =
        BM_ALLOCATE(sizeof(int) * max_efficiencies * n_reactions);
    reactions->efficiencies =
        BM_ALLOCATE(sizeof(double) * max_efficiencies * n_reactions);

    set_value_int_n(0, n_reactions, reactions->n_reactants);
    set_value_int_n(0, n_reactions, reactions->n_products);
    set_value_int_n(0, n_reactions, reactions->n_efficiencies);

    reactions->n_reactions_three = 0;
    reactions->idx_reactions_three = NULL;
    reactions->n_reactions_low = 0;
    reactions->idx_reactions_low = NULL;
    reactions->n_reactions_high = 0;
    reactions->idx_reactions_high = NULL;
    reactions->n_reactions_troe = 0;
    reactions->idx_reactions_troe = NULL;
    reactions->n_reactions_rev = 0;
    reactions->idx_reactions_rev = NULL;
    reactions->n_reactions_rev_arr = 0;
    reactions->idx_reactions_rev_arr = NULL;

    reactions->q = BM_ALLOCATE(sizeof(double) * KINDIM * n_reactions);
    reactions->k = BM_ALLOCATE(sizeof(double) * KINDIM * n_reactions);
    reactions->pr = BM_ALLOCATE(sizeof(double) * n_reactions);

    return reactions;
}

/*******************************************************************************
 * @brief Allocate specii structure
 * @param chemistry
 * @param n_elements
 * @param max_name_length
 * @return elements_t*
 ******************************************************************************/
specii_t *allocate_specii(chemistry_t *chemistry, int n_specii,
                          int max_name_length, int max_reac_points)
{
    chemistry->specii = BM_ALLOCATE(sizeof(specii_t));
    specii_t *specii = chemistry->specii;

    elements_t *elements = chemistry->elements;
    int n_elements = elements->n_elements;

    specii->n_specii = n_specii;
    specii->max_reac_points = max_reac_points;

    specii->symbol =
        allocate_hdf5_string_buffer(n_specii, max_name_length, NULL);
    specii->is_inert = BM_ALLOCATE(sizeof(int) * n_specii);

    specii->composition = BM_ALLOCATE(sizeof(int) * n_elements * n_specii);
    specii->phase = BM_ALLOCATE(sizeof(int) * n_specii);
    specii->bounds = BM_ALLOCATE(sizeof(double) * BOUNDS * n_specii);
    specii->coeff_high = BM_ALLOCATE(sizeof(double) * NASA * n_specii);
    specii->coeff_low = BM_ALLOCATE(sizeof(double) * NASA * n_specii);

    specii->geom = BM_ALLOCATE(sizeof(int) * n_specii);
    specii->pot_lj = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->col_lj = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->dip_mo = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->pol = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->rot_rel = BM_ALLOCATE(sizeof(double) * n_specii);

    specii->n_reac_points = BM_ALLOCATE(sizeof(int) * n_specii);
    specii->reac_points = BM_ALLOCATE(sizeof(int) * max_reac_points * n_specii);
    specii->nu_reac_points =
        BM_ALLOCATE(sizeof(double) * max_reac_points * n_specii);

    set_value_int_n(0, n_specii, specii->n_reac_points);

    specii->molar_mass = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->molecule_weight = BM_ALLOCATE(sizeof(double) * n_specii);
    specii->Rsp = BM_ALLOCATE(sizeof(double) * n_specii);

    specii->omega = BM_ALLOCATE(sizeof(double) * n_specii);

    return specii;
}

/*******************************************************************************
 * @brief Complete chemistry structure
 * @param chemistry
 ******************************************************************************/
void complete_chemistry(chemistry_t *chemistry)
{
    if (chemistry == NULL)
        return;

    elements_t *elements = chemistry->elements;
    int n_elements = elements->n_elements;

    specii_t *specii = chemistry->specii;
    int n_specii = specii->n_specii;

    reactions_t *reactions = chemistry->reactions;
    int n_reactions = reactions->n_reactions;

    for (int i = 0; i < n_specii; ++i)
    {
        specii->molar_mass[i] = 0.0;
        for (int j = 0; j < n_elements; ++j)
            specii->molar_mass[i] +=
                specii->composition[i * n_elements + j] * elements->mass[j];

        specii->molecule_weight[i] = specii->molar_mass[i] / NA;
        specii->Rsp[i] = RM / specii->molar_mass[i];
    }

    reactions->n_reactions_three = 0;
    for (int i = 0; i < n_reactions; ++i)
        if (reactions->type[i] == ReactionThreeBody)
        {
            reactions->n_reactions_three += 1;

            reactions->idx_reactions_three = BM_REALLOCATE(
                reactions->idx_reactions_three,
                sizeof(int) * reactions->n_reactions_three);

            reactions->idx_reactions_three[reactions->n_reactions_three - 1] =
                i;
        }

    reactions->n_reactions_low = 0;
    for (int i = 0; i < n_reactions; ++i)
        if ((reactions->type[i] == ReactionPressure) &&
            (reactions->has_low_arr[i] == 1))
        {
            reactions->n_reactions_low += 1;

            reactions->idx_reactions_low = BM_REALLOCATE(
                reactions->idx_reactions_low,
                sizeof(int) * reactions->n_reactions_low);

            reactions->idx_reactions_low[reactions->n_reactions_low - 1] = i;
        }

    reactions->n_reactions_high = 0;
    for (int i = 0; i < n_reactions; ++i)
        if ((reactions->type[i] == ReactionPressure) &&
            (reactions->has_high_arr[i] == 1))
        {
            reactions->n_reactions_high += 1;

            reactions->idx_reactions_high = BM_REALLOCATE(
                reactions->idx_reactions_high,
                sizeof(int) * reactions->n_reactions_high);

            reactions->idx_reactions_high[reactions->n_reactions_high - 1] = i;
        }

    reactions->n_reactions_troe = 0;
    for (int i = 0; i < n_reactions; ++i)
        if ((reactions->type[i] == ReactionPressure) &&
            (reactions->has_troe[i] == 1))
        {
            reactions->n_reactions_troe += 1;

            reactions->idx_reactions_troe = BM_REALLOCATE(
                reactions->idx_reactions_troe,
                sizeof(int) * reactions->n_reactions_troe);

            reactions->idx_reactions_troe[reactions->n_reactions_troe - 1] = i;
        }

    reactions->n_reactions_rev = 0;
    for (int i = 0; i < n_reactions; ++i)
        if ((reactions->is_reversible[i] == 1) &&
            (reactions->has_rev_arr[i] == 0))
        {
            reactions->n_reactions_rev += 1;

            reactions->idx_reactions_rev = BM_REALLOCATE(
                reactions->idx_reactions_rev,
                sizeof(int) * reactions->n_reactions_rev);

            reactions->idx_reactions_rev[reactions->n_reactions_rev - 1] = i;
        }

    reactions->n_reactions_rev_arr = 0;
    for (int i = 0; i < n_reactions; ++i)
        if ((reactions->is_reversible[i] == 1) &&
            (reactions->has_rev_arr[i] == 1))
        {
            reactions->n_reactions_rev_arr += 1;

            reactions->idx_reactions_rev_arr = BM_REALLOCATE(
                reactions->idx_reactions_rev_arr,
                sizeof(int) * reactions->n_reactions_rev_arr);

            reactions->idx_reactions_rev_arr
                [reactions->n_reactions_rev_arr - 1] = i;
        }
}

/*******************************************************************************
 * @brief Deallocate chemistry structure
 * @param chemistry
 ******************************************************************************/
void deallocate_chemistry(chemistry_t *chemistry)
{
    if (chemistry == NULL)
        return;

    deallocate_elements(chemistry->elements);
    BM_DEALLOCATE(chemistry->elements);

    deallocate_specii(chemistry->specii);
    BM_DEALLOCATE(chemistry->specii);

    deallocate_reactions(chemistry->reactions);
    BM_DEALLOCATE(chemistry->reactions);
}

/*******************************************************************************
 * @brief Deallocate elements structure
 * @param elements
 ******************************************************************************/
void deallocate_elements(elements_t *elements)
{
    if (elements == NULL)
        return;

    deallocate_hdf5_string_buffer(elements->symbol);
    BM_DEALLOCATE(elements->symbol);

    BM_DEALLOCATE(elements->mass);
}

/*******************************************************************************
 * @brief Deallocate reactions structure
 * @param reactions
 ******************************************************************************/
void deallocate_reactions(reactions_t *reactions)
{
    if (reactions == NULL)
        return;

    BM_DEALLOCATE(reactions->type);
    BM_DEALLOCATE(reactions->falloff_species);

    BM_DEALLOCATE(reactions->arr_coeff);
    BM_DEALLOCATE(reactions->is_reversible);

    BM_DEALLOCATE(reactions->n_reactants);
    BM_DEALLOCATE(reactions->reactants);
    BM_DEALLOCATE(reactions->nu_reactants);
    BM_DEALLOCATE(reactions->ord_reactants);

    BM_DEALLOCATE(reactions->n_products);
    BM_DEALLOCATE(reactions->products);
    BM_DEALLOCATE(reactions->nu_products);
    BM_DEALLOCATE(reactions->ord_products);

    BM_DEALLOCATE(reactions->sum_nu);

    BM_DEALLOCATE(reactions->has_rev_arr);
    BM_DEALLOCATE(reactions->rev_arr_coeff);

    BM_DEALLOCATE(reactions->has_high_arr);
    BM_DEALLOCATE(reactions->has_low_arr);
    BM_DEALLOCATE(reactions->adv_arr_coeff);

    BM_DEALLOCATE(reactions->has_troe);
    BM_DEALLOCATE(reactions->n_troe_coeff);
    BM_DEALLOCATE(reactions->troe_coeff);

    BM_DEALLOCATE(reactions->has_efficiencies);
    BM_DEALLOCATE(reactions->n_efficiencies);
    BM_DEALLOCATE(reactions->sp_efficiencies);
    BM_DEALLOCATE(reactions->efficiencies);

    BM_DEALLOCATE(reactions->idx_reactions_three);
    BM_DEALLOCATE(reactions->idx_reactions_low);
    BM_DEALLOCATE(reactions->idx_reactions_high);
    BM_DEALLOCATE(reactions->idx_reactions_troe);
    BM_DEALLOCATE(reactions->idx_reactions_rev);
    BM_DEALLOCATE(reactions->idx_reactions_rev_arr);

    BM_DEALLOCATE(reactions->q);
    BM_DEALLOCATE(reactions->k);
    BM_DEALLOCATE(reactions->pr);
}

/*******************************************************************************
 * @brief Deallocate specii structure
 * @param specii
 ******************************************************************************/
void deallocate_specii(specii_t *specii)
{
    if (specii == NULL)
        return;

    deallocate_hdf5_string_buffer(specii->symbol);
    BM_DEALLOCATE(specii->symbol);

    BM_DEALLOCATE(specii->is_inert);

    BM_DEALLOCATE(specii->composition);
    BM_DEALLOCATE(specii->phase);
    BM_DEALLOCATE(specii->bounds);
    BM_DEALLOCATE(specii->coeff_high);
    BM_DEALLOCATE(specii->coeff_low);

    BM_DEALLOCATE(specii->geom);
    BM_DEALLOCATE(specii->pot_lj);
    BM_DEALLOCATE(specii->col_lj);
    BM_DEALLOCATE(specii->dip_mo);
    BM_DEALLOCATE(specii->pol);
    BM_DEALLOCATE(specii->rot_rel);

    BM_DEALLOCATE(specii->n_reac_points);
    BM_DEALLOCATE(specii->reac_points);
    BM_DEALLOCATE(specii->nu_reac_points);

    BM_DEALLOCATE(specii->molar_mass);
    BM_DEALLOCATE(specii->molecule_weight);
    BM_DEALLOCATE(specii->Rsp);

    BM_DEALLOCATE(specii->omega);
}

/*******************************************************************************
 * @brief Return the species for a given string
 * @param chemistry
 * @param symbol
 * @return int
 ******************************************************************************/
int get_species_index(chemistry_t *chemistry, cstring_t symbol)
{
    specii_t *specii = chemistry->specii;
    int n_specii = specii->n_specii;

    string_t *symbols = specii->symbol;

    for (int i = 0; i < n_specii; ++i)
        if (is_equal(symbol, symbols[i]))
            return i;

    return -1;
}

/*******************************************************************************
 * @brief Print chemistry structure
 * @param chemistry
 ******************************************************************************/
void print_chemistry_info(chemistry_t *chemistry)
{
    if (chemistry == NULL)
        return;

    print_elements(chemistry->elements);
    print_specii(chemistry->specii);
    print_reactions(chemistry->reactions);
}

/*******************************************************************************
 * @brief Print elements structure
 * @param chemistry
 ******************************************************************************/
void print_elements(elements_t *elements)
{
    if (elements == NULL)
        return;

    BM_PRINT("\n");
    printf_r_sep_title(PCSS, "Elements");

    BM_PRINT(PCIN "Number of elements = %d\n", elements->n_elements);

    printf_r_sep(PCSS);
}

/*******************************************************************************
 * @brief Print reactions structure
 * @param chemistry
 ******************************************************************************/
void print_reactions(reactions_t *reactions)
{
    if (reactions == NULL)
        return;

    BM_PRINT("\n");
    printf_r_sep_title(PCSS, "Reactions");

    BM_PRINT(PCIN "Number of reactions = %d\n",
             reactions->n_reactions);
    BM_PRINT(PCIN "Number of three-body reactions = %d\n",
             reactions->n_reactions_three);
    BM_PRINT(PCIN "Number of pressure-dep. reactions (LOW) = %d\n",
             reactions->n_reactions_low);
    BM_PRINT(PCIN "Number of pressure-dep. reactions (HIGH) = %d\n",
             reactions->n_reactions_high);
    BM_PRINT(PCIN "Number of pressure-dep. reactions (TROE) = %d\n",
             reactions->n_reactions_troe);
    BM_PRINT(PCIN "Number of reversible reactions = %d\n",
             reactions->n_reactions_rev);
    BM_PRINT(PCIN "Number of reversible reactions (REV) = %d\n",
             reactions->n_reactions_rev_arr);

    printf_r_sep(PCSS);
}

/*******************************************************************************
 * @brief Print specii structure
 * @param chemistry
 ******************************************************************************/
void print_specii(specii_t *specii)
{
    if (specii == NULL)
        return;

    BM_PRINT("\n");
    printf_r_sep_title(PCSS, "Specii");

    BM_PRINT(PCIN "Number of specii = %d\n", specii->n_specii);

    printf_r_sep(PCSS);
}
