/*******************************************************************************
 * @file restart.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
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

    GET_HDF5_ATTRIBUTE(last_id, "iter", HDF5Int, &iter_restart);
    GET_HDF5_ATTRIBUTE(last_id, "t", HDF5Double, &t_restart);

    phi_restart = ALLOCATE(sizeof(double) * n_variables);
    phi_dt_restart = ALLOCATE(sizeof(double) * n_variables);

    {
        hsize_t dims[1] = {n_variables};
        GET_HDF5_ATTRIBUTE_N(last_id, "phi", HDF5Double, dims[0],
                             phi_restart);
    }

    {
        hsize_t dims[1] = {n_variables};
        GET_HDF5_ATTRIBUTE_N(last_id, "phi_dt", HDF5Double, dims[0],
                             phi_dt_restart);
    }

    if (n_bdf_stages > 0)
    {
        GET_HDF5_ATTRIBUTE(last_id, "n_stages", HDF5Int, &n_stages_restart);
        CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

        phi_old_restart = ALLOCATE(sizeof(double *) * n_stages_restart);
        for (int i = 0; i < n_stages_restart; ++i)
            phi_old_restart[i] = ALLOCATE(sizeof(double) * n_variables);

        for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
        {
            char iter_string[11];
            sprintf(iter_string, "%d", i_stage);
            string_t tmp = allocate_strcat("phi_old:", iter_string);

            hsize_t dims[1] = {n_variables};
            GET_HDF5_DATASET_N(last_id, tmp, HDF5Double, dims[0],
                               phi_old_restart[i_stage]);

            DEALLOCATE(tmp);
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

    DEALLOCATE(phi_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define()
{
    REGISTER_INITIALIZE_ROUTINE(restart_initialize);
    REGISTER_FINALIZE_ROUTINE(restart_finalize);

    SET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart,
                  "The flag to start from restart", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize()
{
    DEALLOCATE(phi_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize()
{
    GET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart);

    if (use_restart == 1)
        read_restart_data();
    else
        create_file_header();
}
