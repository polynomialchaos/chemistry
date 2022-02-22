/*******************************************************************************
 * @file restart_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef RESTART_MODULE_H
#define RESTART_MODULE_H

#include "chemistry/reactor0d_module.h"

extern int use_restart;

extern int iter_restart;
extern double t_restart;

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define();

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize();

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize();

#endif /* RESTART_MODULE_H */