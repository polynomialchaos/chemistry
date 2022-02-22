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
    chemistry_t *tmp = ALLOCATE(sizeof(chemistry_t));

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
    chemistry->elements = ALLOCATE(sizeof(elements_t));
    elements_t *elements = chemistry->elements;

    elements->n_elements = n_elements;

    elements->symbol =
        allocate_hdf5_string_buffer(n_elements, max_name_length, NULL);
    elements->mass = ALLOCATE(sizeof(double) * n_elements);

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
    chemistry->reactions = ALLOCATE(sizeof(reactions_t));
    reactions_t *reactions = chemistry->reactions;

    reactions->n_reactions = n_reactions;
    reactions->max_reactants = max_reactants;
    reactions->max_products = max_products;
    reactions->max_troe_coeff = max_troe_coeff;
    reactions->max_efficiencies = max_efficiencies;

    reactions->type = ALLOCATE(sizeof(int) * n_reactions);
    reactions->falloff_species = ALLOCATE(sizeof(int) * n_reactions);

    reactions->arr_coeff = ALLOCATE(sizeof(double) * ARR * n_reactions);
    reactions->is_reversible = ALLOCATE(sizeof(int) * n_reactions);

    reactions->n_reactants = ALLOCATE(sizeof(int) * n_reactions);
    reactions->reactants = ALLOCATE(sizeof(int) * max_reactants * n_reactions);
    reactions->nu_reactants =
        ALLOCATE(sizeof(double) * max_reactants * n_reactions);
    reactions->ord_reactants =
        ALLOCATE(sizeof(double) * max_reactants * n_reactions);

    reactions->n_products = ALLOCATE(sizeof(int) * n_reactions);
    reactions->products = ALLOCATE(sizeof(int) * max_products * n_reactions);
    reactions->nu_products =
        ALLOCATE(sizeof(double) * max_products * n_reactions);
    reactions->ord_products =
        ALLOCATE(sizeof(double) * max_products * n_reactions);

    reactions->sum_nu = ALLOCATE(sizeof(double) * n_reactions);

    reactions->has_rev_arr = ALLOCATE(sizeof(int) * n_reactions);
    reactions->rev_arr_coeff = ALLOCATE(sizeof(double) * ARR * n_reactions);

    reactions->has_high_arr = ALLOCATE(sizeof(int) * n_reactions);
    reactions->has_low_arr = ALLOCATE(sizeof(int) * n_reactions);
    reactions->adv_arr_coeff = ALLOCATE(sizeof(double) * ARR * n_reactions);

    reactions->has_troe = ALLOCATE(sizeof(int) * n_reactions);
    reactions->n_troe_coeff = ALLOCATE(sizeof(int) * n_reactions);
    reactions->troe_coeff =
        ALLOCATE(sizeof(double) * max_troe_coeff * n_reactions);

    reactions->has_efficiencies = ALLOCATE(sizeof(int) * n_reactions);
    reactions->n_efficiencies = ALLOCATE(sizeof(int) * n_reactions);
    reactions->sp_efficiencies =
        ALLOCATE(sizeof(int) * max_efficiencies * n_reactions);
    reactions->efficiencies =
        ALLOCATE(sizeof(double) * max_efficiencies * n_reactions);

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

    reactions->q = ALLOCATE(sizeof(double) * KINDIM * n_reactions);
    reactions->k = ALLOCATE(sizeof(double) * KINDIM * n_reactions);
    reactions->pr = ALLOCATE(sizeof(double) * n_reactions);

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
    chemistry->specii = ALLOCATE(sizeof(specii_t));
    specii_t *specii = chemistry->specii;

    elements_t *elements = chemistry->elements;
    int n_elements = elements->n_elements;

    specii->n_specii = n_specii;
    specii->max_reac_points = max_reac_points;

    specii->symbol =
        allocate_hdf5_string_buffer(n_specii, max_name_length, NULL);
    specii->is_inert = ALLOCATE(sizeof(int) * n_specii);

    specii->composition = ALLOCATE(sizeof(int) * n_elements * n_specii);
    specii->phase = ALLOCATE(sizeof(int) * n_specii);
    specii->bounds = ALLOCATE(sizeof(double) * BOUNDS * n_specii);
    specii->coeff_high = ALLOCATE(sizeof(double) * NASA * n_specii);
    specii->coeff_low = ALLOCATE(sizeof(double) * NASA * n_specii);

    specii->geom = ALLOCATE(sizeof(int) * n_specii);
    specii->pot_lj = ALLOCATE(sizeof(double) * n_specii);
    specii->col_lj = ALLOCATE(sizeof(double) * n_specii);
    specii->dip_mo = ALLOCATE(sizeof(double) * n_specii);
    specii->pol = ALLOCATE(sizeof(double) * n_specii);
    specii->rot_rel = ALLOCATE(sizeof(double) * n_specii);

    specii->n_reac_points = ALLOCATE(sizeof(int) * n_specii);
    specii->reac_points = ALLOCATE(sizeof(int) * max_reac_points * n_specii);
    specii->nu_reac_points =
        ALLOCATE(sizeof(double) * max_reac_points * n_specii);

    set_value_int_n(0, n_specii, specii->n_reac_points);

    specii->molar_mass = ALLOCATE(sizeof(double) * n_specii);
    specii->molecule_weight = ALLOCATE(sizeof(double) * n_specii);
    specii->Rsp = ALLOCATE(sizeof(double) * n_specii);

    specii->omega = ALLOCATE(sizeof(double) * n_specii);

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

            reactions->idx_reactions_three = REALLOCATE(
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

            reactions->idx_reactions_low = REALLOCATE(
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

            reactions->idx_reactions_high = REALLOCATE(
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

            reactions->idx_reactions_troe = REALLOCATE(
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

            reactions->idx_reactions_rev = REALLOCATE(
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

            reactions->idx_reactions_rev_arr = REALLOCATE(
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
    DEALLOCATE(chemistry->elements);

    deallocate_specii(chemistry->specii);
    DEALLOCATE(chemistry->specii);

    deallocate_reactions(chemistry->reactions);
    DEALLOCATE(chemistry->reactions);
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
    DEALLOCATE(elements->symbol);

    DEALLOCATE(elements->mass);
}

/*******************************************************************************
 * @brief Deallocate reactions structure
 * @param reactions
 ******************************************************************************/
void deallocate_reactions(reactions_t *reactions)
{
    if (reactions == NULL)
        return;

    DEALLOCATE(reactions->type);
    DEALLOCATE(reactions->falloff_species);

    DEALLOCATE(reactions->arr_coeff);
    DEALLOCATE(reactions->is_reversible);

    DEALLOCATE(reactions->n_reactants);
    DEALLOCATE(reactions->reactants);
    DEALLOCATE(reactions->nu_reactants);
    DEALLOCATE(reactions->ord_reactants);

    DEALLOCATE(reactions->n_products);
    DEALLOCATE(reactions->products);
    DEALLOCATE(reactions->nu_products);
    DEALLOCATE(reactions->ord_products);

    DEALLOCATE(reactions->sum_nu);

    DEALLOCATE(reactions->has_rev_arr);
    DEALLOCATE(reactions->rev_arr_coeff);

    DEALLOCATE(reactions->has_high_arr);
    DEALLOCATE(reactions->has_low_arr);
    DEALLOCATE(reactions->adv_arr_coeff);

    DEALLOCATE(reactions->has_troe);
    DEALLOCATE(reactions->n_troe_coeff);
    DEALLOCATE(reactions->troe_coeff);

    DEALLOCATE(reactions->has_efficiencies);
    DEALLOCATE(reactions->n_efficiencies);
    DEALLOCATE(reactions->sp_efficiencies);
    DEALLOCATE(reactions->efficiencies);

    DEALLOCATE(reactions->idx_reactions_three);
    DEALLOCATE(reactions->idx_reactions_low);
    DEALLOCATE(reactions->idx_reactions_high);
    DEALLOCATE(reactions->idx_reactions_troe);
    DEALLOCATE(reactions->idx_reactions_rev);
    DEALLOCATE(reactions->idx_reactions_rev_arr);

    DEALLOCATE(reactions->q);
    DEALLOCATE(reactions->k);
    DEALLOCATE(reactions->pr);
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
    DEALLOCATE(specii->symbol);

    DEALLOCATE(specii->is_inert);

    DEALLOCATE(specii->composition);
    DEALLOCATE(specii->phase);
    DEALLOCATE(specii->bounds);
    DEALLOCATE(specii->coeff_high);
    DEALLOCATE(specii->coeff_low);

    DEALLOCATE(specii->geom);
    DEALLOCATE(specii->pot_lj);
    DEALLOCATE(specii->col_lj);
    DEALLOCATE(specii->dip_mo);
    DEALLOCATE(specii->pol);
    DEALLOCATE(specii->rot_rel);

    DEALLOCATE(specii->n_reac_points);
    DEALLOCATE(specii->reac_points);
    DEALLOCATE(specii->nu_reac_points);

    DEALLOCATE(specii->molar_mass);
    DEALLOCATE(specii->molecule_weight);
    DEALLOCATE(specii->Rsp);

    DEALLOCATE(specii->omega);
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

    PRINTF("\n");
    printf_r_sep_title(PCSS, "Elements");

    PRINTF(PCIN "Number of elements = %d\n", elements->n_elements);

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

    PRINTF("\n");
    printf_r_sep_title(PCSS, "Reactions");

    PRINTF(PCIN "Number of reactions = %d\n",
        reactions->n_reactions);
    PRINTF(PCIN "Number of three-body reactions = %d\n",
        reactions->n_reactions_three);
    PRINTF(PCIN "Number of pressure-dep. reactions (LOW) = %d\n",
        reactions->n_reactions_low);
    PRINTF(PCIN "Number of pressure-dep. reactions (HIGH) = %d\n",
        reactions->n_reactions_high);
    PRINTF(PCIN "Number of pressure-dep. reactions (TROE) = %d\n",
        reactions->n_reactions_troe);
    PRINTF(PCIN "Number of reversible reactions = %d\n",
        reactions->n_reactions_rev);
    PRINTF(PCIN "Number of reversible reactions (REV) = %d\n",
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

    PRINTF("\n");
    printf_r_sep_title(PCSS, "Specii");

    PRINTF(PCIN "Number of specii = %d\n", specii->n_specii);

    printf_r_sep(PCSS);
}
