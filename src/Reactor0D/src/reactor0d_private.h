/*******************************************************************************
 * @file reactor0d_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef REACTOR0D_PRIVATE_H
#define REACTOR0D_PRIVATE_H

#include "chemistry/reactor0d_module.h"

extern string_t title;

/*******************************************************************************
 * @brief Define reactor0d
 ******************************************************************************/
void reactor0d_define();

/*******************************************************************************
 * @brief Finalize reactor0d
 ******************************************************************************/
void reactor0d_finalize();

/*******************************************************************************
 * @brief Initialize reactor0d
 ******************************************************************************/
void reactor0d_initialize();

#endif /* REACTOR0D_PRIVATE_H */