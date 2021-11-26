/*******************************************************************************
 * @file output_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include "chemistry/reactor0d_module.h"

extern int do_output_data;
extern int i_output_data;
extern string_t output_file;

/*******************************************************************************
 * @brief Create a file header
 ******************************************************************************/
void create_file_header();

/*******************************************************************************
 * @brief Define output
 ******************************************************************************/
void output_define();

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void output_finalize();

/*******************************************************************************
 * @brief Initialize output
 ******************************************************************************/
void output_initialize();

/*******************************************************************************
 * @brief Write data to file
 * @param iter
 * @param t
 ******************************************************************************/
void write_output(int iter, double t);

#endif /* OUTPUT_MODULE_H */