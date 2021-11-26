/*******************************************************************************
 * @file implicit_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef IMPLICIT_MODULE_H
#define IMPLICIT_MODULE_H

#include "chemistry/reactor0d_module.h"

extern int implicit_active;
extern int n_iter_inner;
extern int n_iter_lsoe;
extern int n_bdf_stages;
extern double **phi_old;

/*******************************************************************************
 * @brief Numerical jacobian routine
 * @param n_var
 * @param n_cells
 ******************************************************************************/
void calc_jacobian_numerical(int n_var);

/*******************************************************************************
 * @brief Define implicit timedisc
 ******************************************************************************/
void implicit_define();

/*******************************************************************************
 * @brief Finalize implicit timedisc
 ******************************************************************************/
void implicit_finalize();

/*******************************************************************************
 * @brief Initialize implicit timedisc
 ******************************************************************************/
void implicit_initialize();

/*******************************************************************************
 * @brief Implicit time discretizazion routine (Newton)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_newton(int iter, double t, double dt);

/*******************************************************************************
 * @brief Matrix vector routine (called by solver)
 * @param x
 * @param b
 * @param n_var
 * @param m
 * @return int
 ******************************************************************************/
int matrix_vector_numerical(double *x, double *b, size_t n_var, size_t m);

#endif /* IMPLICIT_MODULE_H */