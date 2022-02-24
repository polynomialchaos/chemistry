/*******************************************************************************
 * @file chemistry.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include <stdio.h>
#include "chemistry/chemistry_module.h"

/*******************************************************************************
 * @brief Main function
 * @param argc
 * @param argv
 * @return int
 ******************************************************************************/
int main(int argc, char **argv)
{
    global_initialize(argc, argv, BC_FALSE, BC_FALSE, BC_FALSE, BC_FALSE);

    string_t chem_file = argc > 1 ? argv[1] : "../gri3.mech.h5";
    chemistry_t *chemistry = read_chemistry_data(chem_file);

    print_chemistry_info(chemistry);
    deallocate_chemistry(chemistry);
    BM_DEALLOCATE(chemistry);

    check_abort(1);
    return 0;
}