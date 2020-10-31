//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include "reactor0d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
extern int do_output_data;
extern int i_output_data;
extern string_t output_file;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_define();

void create_file_header();
void write_output( int iter, double t );

#endif /* OUTPUT_MODULE_H */