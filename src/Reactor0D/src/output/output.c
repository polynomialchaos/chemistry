/*******************************************************************************
 * @file output.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <string.h>
#include "reactor0d_private.h"
#include "reactor/reactor_module.h"
#include "output_module.h"
#include "timedisc/implicit/implicit_module.h"

int i_output_data = -1;
int do_output_data = 0;
string_t output_file = NULL;

/*******************************************************************************
 * @brief Create a file header
 ******************************************************************************/
void create_file_header()
{
    hid_t file_id = create_hdf5_file(output_file);

    SET_HDF5_ATTRIBUTE(file_id, "n_variables", HDF5Int, &n_variables);

    {
        size_t max_len = strlen_n(variables, n_variables) + 1;
        string_t *tmp =
            allocate_hdf5_string_buffer(n_variables, max_len, variables);

        hsize_t dims[1] = {n_variables};
        SET_HDF5_DATASET_N(file_id, "variables", HDF5String, tmp, dims[0]);

        deallocate_hdf5_string_buffer(tmp);
        DEALLOCATE(tmp);
    }

    hid_t group_id = create_hdf5_group(file_id, "SOLUTIONS");
    close_hdf5_group(group_id);

    close_hdf5_file(file_id);
}

/*******************************************************************************
 * @brief Write data to file
 * @param iter
 * @param t
 ******************************************************************************/
void write_output(int iter, double t)
{
    if (is_valid_hdf5_file(output_file) == 0)
        create_file_header();

    char iter_string[256];
    sprintf(iter_string, "%09d", iter);

    hid_t file_id = open_hdf5_file(output_file);

    hid_t group_id = open_hdf5_group(file_id, "SOLUTIONS");

    hid_t solution_id = create_hdf5_group(group_id, iter_string);

    SET_HDF5_ATTRIBUTE(solution_id, "iter", HDF5Int, &iter);
    SET_HDF5_ATTRIBUTE(solution_id, "t", HDF5Double, &t);

    {
        hsize_t dims[1] = {n_variables};
        SET_HDF5_DATASET_N(solution_id, "phi", HDF5Double, phi, dims[0]);
    }

    {
        hsize_t dims[1] = {n_variables};
        SET_HDF5_DATASET_N(solution_id, "phi_dt", HDF5Double, phi_dt, dims[0]);
    }

    if (n_bdf_stages > 0)
    {
        SET_HDF5_ATTRIBUTE(solution_id, "n_stages", HDF5Int, &n_bdf_stages);

        for (int i_stage = 0; i_stage < n_bdf_stages; ++i_stage)
        {
            char stage_string[256];
            sprintf(stage_string, "%d", i_stage);
            string_t tmp = allocate_strcat("phi_old:", stage_string);

            hsize_t dims[1] = {n_variables};

            SET_HDF5_DATASET_N(solution_id, tmp, HDF5Double,
                               phi_old[i_stage], dims[0]);

            DEALLOCATE(tmp);
        }
    }

    close_hdf5_group(solution_id);

    close_hdf5_group(group_id);

    if (exists_hdf5_link(file_id, "SOLUTION"))
        delete_hdf5_link(file_id, "SOLUTION");
    string_t tmp = allocate_strcat("SOLUTIONS/", iter_string);
    create_hdf5_soft_link(file_id, "SOLUTION", tmp);
    DEALLOCATE(tmp);

    close_hdf5_file(file_id);
}

/*******************************************************************************
 * @brief Define output
 ******************************************************************************/
void output_define()
{
    REGISTER_INITIALIZE_ROUTINE(output_initialize);
    REGISTER_FINALIZE_ROUTINE(output_finalize);

    SET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data,
                  "The output file frequency " \
                  "(-1 ... first/solutions/last, 0 ... disable)",
                  NULL, 0);
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void output_finalize()
{
    DEALLOCATE(output_file);
}

/*******************************************************************************
 * @brief Initialize output
 ******************************************************************************/
void output_initialize()
{
    GET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data);
    do_output_data = (i_output_data != 0);

    output_file = allocate_strcat(title, ".h5");
}
