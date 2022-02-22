/*******************************************************************************
 * @file reactor_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef REACTOR_MODULE_H
#define REACTOR_MODULE_H

#include "chemistry/reactor0d_module.h"

typedef void (*void_reactor_ft)(double t);
extern void_reactor_ft reactor_function_pointer;

extern const int i_Y0;
extern int n_variables;
extern string_t *variables;

extern double *phi;
extern double *phi_dt;
extern double *phi_bounds;

/*******************************************************************************
 * @brief Return the temperature equation source term
 ******************************************************************************/
double calc_dT_dt(double *dY_dt, double T_old);

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void reactor_define();

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void reactor_finalize();

/*******************************************************************************
 * @brief Initialize reactor
 ******************************************************************************/
void reactor_initialize();

#endif /* REACTOR_MODULE_H */