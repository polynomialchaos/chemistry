//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef RESTART_MODULE_H
#define RESTART_MODULE_H

#include "chemistry/reactor0d_module.h"




extern int use_restart;

extern int iter_restart;
extern double t_restart;


void restart_define();

#endif /* RESTART_MODULE_H */