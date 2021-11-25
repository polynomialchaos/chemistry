/*******************************************************************************
 * @file reactor0d.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "reactor0d_private.h"
#include "reactor/reactor_module.h"
#include "analyze/analyze_module.h"
#include "output/output_module.h"
#include "restart/restart_module.h"
#include "timedisc/timedisc_module.h"

string_t title = NULL;

/*******************************************************************************
 * @brief Main function
 * @param argc
 * @param argv
 * @return int
 ******************************************************************************/
int main(int argc, string_t *argv)
{
    /* define the program structure */
    reactor0d_define();
    reactor_define();
    analyze_define();
    output_define();
    timedisc_define();
    restart_define();

    /* call the global initialize routine */
    global_initialize(argc, argv, BFLS, BFLS, BFLS, BTRU);

    /* calculation */
    PRINTF("\n");
    printf_r_sep_title('=', "Calculation");

    timedisc();

    printf_r_sep('=');

    /* end the program */
    check_abort(1);
    return 1;
}

/*******************************************************************************
 * @brief Define reactor0d
 ******************************************************************************/
void reactor0d_define()
{
    REGISTER_INITIALIZE_ROUTINE(reactor0d_initialize);
    REGISTER_FINALIZE_ROUTINE(reactor0d_finalize);

    string_t tmp = "untitled";
    SET_PARAMETER("General/title", StringParameter, &tmp,
                  "The project title", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize reactor0d
 ******************************************************************************/
void reactor0d_finalize()
{
    DEALLOCATE(title);
}

/*******************************************************************************
 * @brief Initialize reactor0d
 ******************************************************************************/
void reactor0d_initialize()
{
    GET_PARAMETER("General/title", StringParameter, &title);
}