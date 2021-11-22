//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef REACTOR_MODULE_H
#define REACTOR_MODULE_H

#include "chemistry/reactor0d_module.h"




typedef void (*void_reactor_fp_t)(double t);
extern void_reactor_fp_t reactor_function_pointer;

extern const int i_Y0;
extern int n_variables;
extern string_t *variables;

extern double *phi;
extern double *phi_dt;
extern double *phi_bounds;


void reactor_define();

#endif /* REACTOR_MODULE_H */