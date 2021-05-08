//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "reactor0d_private.h"
#include "reactor/reactor_module.h"
#include "output_module.h"
#include "timedisc/implicit/implicit_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int i_output_data = -1;
int do_output_data = 0;
string_t output_file = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_initialize();
void output_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_define()
{
    register_initialize_routine(output_initialize);
    register_finalize_routine(output_finalize);

    set_parameter("Output/i_output_data", ParameterDigit, &i_output_data,
                  "The output file frequency  (-1 ... first/solutions/last, 0 ... disable)", NULL, 0);
}

void output_initialize()
{
    get_parameter("Output/i_output_data", ParameterDigit, &i_output_data);
    do_output_data = (i_output_data != 0);

    output_file = allocate_strcat(title, ".h5");
}

void output_finalize()
{
    deallocate(output_file);
}

void create_file_header()
{
    hid_t file_id = create_hdf5_file(output_file);

    set_hdf5_attribute(file_id, "n_variables", HDF5Int, &n_variables);

    {
        size_t max_len = 0;
        for (int i = 0; i < n_variables; i++)
            max_len = u_max(max_len, strlen(variables[i]));

        string_t *tmp = allocate_hdf5_string_buffer(n_variables, max_len + 1);
        for (int i = 0; i < n_variables; i++)
            strcpy(tmp[i], variables[i]);

        hsize_t dims[2] = {n_variables, max_len + 1};
        set_hdf5_dataset_n(file_id, "variables", HDF5String, tmp[0], dims);

        deallocate_hdf5_string_buffer(&tmp);
    }

    hid_t group_id = create_hdf5_group(file_id, "SOLUTIONS");
    close_hdf5_group(group_id);

    close_hdf5_file(file_id);
}

void write_output(int iter, double t)
{
    if (is_valid_hdf5_file(output_file) == 0)
        create_file_header();

    char iter_string[256];
    sprintf(iter_string, "%09d", iter);

    hid_t file_id = open_hdf5_file(output_file);

    hid_t group_id = open_hdf5_group(file_id, "SOLUTIONS");

    hid_t solution_id = create_hdf5_group(group_id, iter_string);

    set_hdf5_attribute(solution_id, "iter", HDF5Int, &iter);
    set_hdf5_attribute(solution_id, "t", HDF5Double, &t);

    {
        hsize_t dims[1] = {n_variables};
        set_hdf5_dataset_n(solution_id, "phi", HDF5Double, phi, dims);
    }

    {
        hsize_t dims[1] = {n_variables};
        set_hdf5_dataset_n(solution_id, "phi_dt", HDF5Double, phi_dt, dims);
    }

    if (n_bdf_stages > 0)
    {
        set_hdf5_attribute(solution_id, "n_stages", HDF5Int, &n_bdf_stages);

        for (int i_stage = 0; i_stage < n_bdf_stages; i_stage++)
        {
            char stage_string[256];
            sprintf(stage_string, "%d", i_stage);
            string_t tmp = allocate_strcat("phi_old:", stage_string);

            hsize_t dims[1] = {n_variables};

            set_hdf5_dataset_n(solution_id, tmp, HDF5Double, phi_old[i_stage], dims);

            deallocate(tmp);
        }
    }

    close_hdf5_group(solution_id);

    close_hdf5_group(group_id);

    if (exists_hdf5_link(file_id, "SOLUTION"))
        delete_hdf5_link(file_id, "SOLUTION");
    string_t tmp = allocate_strcat("SOLUTIONS/", iter_string);
    create_hdf5_soft_link(file_id, "SOLUTION", tmp);
    deallocate(tmp);

    close_hdf5_file(file_id);
}