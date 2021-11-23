/*******************************************************************************
 * @file chemistry.c
 * @author Florian Eigentler
 * @brief
 * @version 0.1
 * @date 2021-11-08
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <stdio.h>
#include "chemistry/chemistry_module.h"

/*******************************************************************************
 * @brief Main test function
 * @param argc
 * @param argv
 * @return int
 ******************************************************************************/
int main(int argc, char **argv)
{
    global_initialize(argc, argv, BFLS, BFLS, BFLS, BFLS);

    string_t chem_file = argc > 1 ? argv[1] : "../gri3.mech.h5";
    chemistry_t *chemistry = read_chemistry_data(chem_file);

    print_chemistry_info(chemistry);
    deallocate_chemistry(chemistry);
    DEALLOCATE(chemistry);

    check_abort(1);
    return 0;
}