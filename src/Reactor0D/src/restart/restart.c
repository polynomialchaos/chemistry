/*******************************************************************************
 * @file restart.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "restart_module.h"
#include "reactor/reactor_module.h"
#include "output/output_module.h"
#include "timedisc/implicit/implicit_module.h"

int use_restart = 0;

int iter_restart = 0;
double t_restart = 0;

double *phi_restart = NULL;
double *phi_dt_restart = NULL;
double **phi_old_restart = NULL;

int n_stages_restart = 0;

/*******************************************************************************
 * @brief Read restart data
 ******************************************************************************/
void read_restart_data()
{
    hid_t file_id = open_hdf5_file(output_file);

    hid_t last_id = open_hdf5_group(file_id, "SOLUTION");

    BM_GET_HDF5_ATTRIBUTE(last_id, "iter", HDF5Int, &iter_restart);
    BM_GET_HDF5_ATTRIBUTE(last_id, "t", HDF5Double, &t_restart);

    phi_restart = BM_ALLOCATE(sizeof(double) * n_variables);
    phi_dt_restart = BM_ALLOCATE(sizeof(double) * n_variables);

    {
        hsize_t dims[1] = {n_variables};
        BM_GET_HDF5_ATTRIBUTE_N(last_id, "phi", HDF5Double, dims[0],
                                phi_restart);
    }

    {
        hsize_t dims[1] = {n_variables};
        BM_GET_HDF5_ATTRIBUTE_N(last_id, "phi_dt", HDF5Double, dims[0],
                                phi_dt_restart);
    }

    if (n_bdf_stages > 0)
    {
        BM_GET_HDF5_ATTRIBUTE(last_id, "n_stages", HDF5Int, &n_stages_restart);
        BM_CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

        phi_old_restart = BM_ALLOCATE(sizeof(double *) * n_stages_restart);
        for (int i = 0; i < n_stages_restart; ++i)
            phi_old_restart[i] = BM_ALLOCATE(sizeof(double) * n_variables);

        for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
        {
            char iter_string[11];
            sprintf(iter_string, "%d", i_stage);
            string_t tmp = allocate_strcat("phi_old:", iter_string);

            hsize_t dims[1] = {n_variables};
            BM_GET_HDF5_DATASET_N(last_id, tmp, HDF5Double, dims[0],
                                  phi_old_restart[i_stage]);

            BM_DEALLOCATE(tmp);
        }
    }

    close_hdf5_group(last_id);

    close_hdf5_file(file_id);

    copy_n(phi_restart, n_variables, phi);
    copy_n(phi_dt_restart, n_variables, phi_dt);

    for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
    {
        copy_n(phi_old_restart[i_stage], n_variables, phi_old[i_stage]);
    }

    BM_DEALLOCATE(phi_restart);
    BM_DEALLOCATE(phi_dt_restart);
    BM_DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(restart_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(restart_finalize);

    BM_SET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart,
                     "The flag to start from restart", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize()
{
    BM_DEALLOCATE(phi_restart);
    BM_DEALLOCATE(phi_dt_restart);
    BM_DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize()
{
    BM_GET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart);

    if (use_restart == 1)
        read_restart_data();
    else
        create_file_header();
}
