/*******************************************************************************
 * @file explicit_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef EXPLICIT_MODULE_H
#define EXPLICIT_MODULE_H

#include "chemistry/reactor0d_module.h"

extern int explicit_active;

/*******************************************************************************
 * @brief Define explicit timedisc
 ******************************************************************************/
void explicit_define();

/*******************************************************************************
 * @brief Finalize explicit timedisc
 ******************************************************************************/
void explicit_finalize();

/*******************************************************************************
 * @brief Initialize explicit timedisc
 ******************************************************************************/
void explicit_initialize();

/*******************************************************************************
 * @brief Explicit time discretizazion routine (LSERKW2)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_lserkw2(int iter, double t, double dt);

#endif /* EXPLICIT_MODULE_H */