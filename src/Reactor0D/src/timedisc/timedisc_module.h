/*******************************************************************************
 * @file timedisc_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef TIMEDISC_MODULE_H
#define TIMEDISC_MODULE_H

#include "chemistry/reactor0d_module.h"

typedef void (*void_timestep_ft)(int iter, double t, double dt);
extern void_timestep_ft time_step_function_pointer;

typedef double (*double_calc_timestep_ft)();
extern double_calc_timestep_ft calc_time_step_function_pointer;

extern int is_viscous_dt;
extern int is_transient;

/*******************************************************************************
 * @brief Time discretizazion routine
 ******************************************************************************/
void timedisc();

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define();

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void timedisc_finalize();

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize();

#endif /* TIMEDISC_MODULE_H */