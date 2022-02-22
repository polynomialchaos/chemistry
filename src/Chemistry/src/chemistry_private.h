/*******************************************************************************
 * @file chemistry_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef CHEMISTRY_PRIVATE_H
#define CHEMISTRY_PRIVATE_H

#include "chemistry/chemistry_module.h"

/*******************************************************************************
 * @brief Reaction type enumeration
 ******************************************************************************/
typedef enum ReactionType
{
    ReactionDefault,   /** Default reaction */
    ReactionThreeBody, /** Three-body reaction */
    ReactionPressure,  /** Pressure dependent reaction */
    _reaction_type_max
} reaction_type_t;

/*******************************************************************************
 * @brief Allocate chemistry structure
 * @return chemistry_t*
 ******************************************************************************/
chemistry_t *allocate_chemistry();

/*******************************************************************************
 * @brief Allocate elements structure
 * @param chemistry
 * @param n_elements
 * @param max_name_length
 * @return elements_t*
 ******************************************************************************/
elements_t *allocate_elements(chemistry_t *chemistry,
                              int n_elements, int max_name_length);

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
                                int max_troe_coeff, int max_efficiencies);

/*******************************************************************************
 * @brief Allocate specii structure
 * @param chemistry
 * @param n_elements
 * @param max_name_length
 * @return elements_t*
 ******************************************************************************/
specii_t *allocate_specii(chemistry_t *chemistry, int n_specii,
                          int max_name_length, int max_reac_points);

/*******************************************************************************
 * @brief Calculate the Arrhenius function
 * @param T
 * @param pr
 * @param coeff
 * @return double
 ******************************************************************************/
double calc_troe(double T, double pr, double *coeff);

/*******************************************************************************
 * @brief Complete chemistry structure
 * @param chemistry
 ******************************************************************************/
void complete_chemistry(chemistry_t *chemistry);

/*******************************************************************************
 * @brief Deallocate chemistry structure
 * @param elements
 ******************************************************************************/
void deallocate_elements(elements_t *elements);

/*******************************************************************************
 * @brief Deallocate chemistry structure
 * @param reactions
 ******************************************************************************/
void deallocate_reactions(reactions_t *reactions);

/*******************************************************************************
 * @brief Deallocate chemistry structure
 * @param specii
 ******************************************************************************/
void deallocate_specii(specii_t *specii);

/*******************************************************************************
 * @brief Print elements structure
 * @param chemistry
 ******************************************************************************/
void print_elements(elements_t *elements);

/*******************************************************************************
 * @brief Print reactions structure
 * @param chemistry
 ******************************************************************************/
void print_reactions(reactions_t *reactions);

/*******************************************************************************
 * @brief Print specii structure
 * @param chemistry
 ******************************************************************************/
void print_specii(specii_t *specii);

#endif /* CHEMISTRY_PRIVATE_H */