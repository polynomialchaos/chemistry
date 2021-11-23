/*******************************************************************************
 * @file xxx.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef RESTART_MODULE_H
#define RESTART_MODULE_H

#include "chemistry/reactor0d_module.h"




extern int use_restart;

extern int iter_restart;
extern double t_restart;


void restart_define();

#endif /* RESTART_MODULE_H */