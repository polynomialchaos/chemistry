//##################################################################################################################################
// Chemistry - Finite rate chemistry library and solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef CHEMISTRY_MODULE_H
#define CHEMISTRY_MODULE_H

#include "c_module.h"
#include "utilities_module.h"
#include "parameter_module.h"
#include "math_module.h"
#include "mpi_module.h"
#include "global_module.h"
#include "hdf5_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------
#define NA 6.022140857e23   // Avogadro number, mol-1
#define KB 1.38064852e-23   // Boltzmann constant, m2 kg s-2 K-1
#define RM NA * KB          // universal gas constant, mol-1 m2 kg s-2 K-1
#define P0 101325.0         // standard pressure, Pa
#define T0 273.15           // standard temperature, K

#define BOUNDS 3            // number of temperature bounds
#define NASA 7              // number of NASA polynomial coefficients
#define ARR 3               // number of Arrhenius coefficients
#define TROE 5              // number of Troe pressure function coefficients
#define KINDIM 2            // nmuber of kinetic dimensions (forward/backward)

#define YONE 1e-6           // maximum tolerance for sum y - 1
#define YSMALL 1e-100       // minimum species concentration

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
typedef struct Elements
{
    int          n_elements;            // number of elements

    string_t    *symbol;                // symbol
    double      *mass;                  // mass, kg mol-1
} Elements_t;

typedef struct Specii
{
    int          n_specii;              // number of specii
    int          max_reac_points;       // maximum number of reaction points

    string_t    *symbol;                // symbol
    int         *is_inert;              // flag for inert species

    int         *composition;           // elementary composition
    int         *phase;                 // phase
    double      *bounds;                // temperature bounds for NASA polynomials
    double      *coeff_high;            // coefficients for high temperature range
    double      *coeff_low;             // coefficients for low temperature range

    int         *geom;                  // geometrical configuration
    double      *pot_lj;                // lennard-Jones potential
    double      *col_lj;                // lennard-Jones collision diameter
    double      *dip_mo;                // dipole momentum
    double      *pol;                   // polarizeability
    double      *rot_rel;               // rotational relaxation

    int         *n_reac_points;         // number of reaction points
    int         *reac_points;           // reaction points
    double      *nu_reac_points;        // reaction points nu

    double      *molar_mass;            // molar mass of species, kg mol-1
    double      *molecule_weight;       // molecular mass of species, kg
    double      *Rsp;                   // specific gas constant, J kg-1 K-1

    double      *omega;                 // production rates (kg m-3 s-1)
} Specii_t;

typedef struct Reactions
{
    int          n_reactions;           // number of reactions
    int          max_reactants;         // maximum number of reactants
    int          max_products;          // maximum number of products
    int          max_troe_coeff;        // maximum number of troe coefficients
    int          max_efficiencies;      // maximum number of efficiencies

    int         *type;                  // reaction type
    int         *falloff_species;       // fall off species

    double      *arr_coeff;             // Arrhenius coefficient
    int         *is_reversible;         // flag for reversibility

    int         *n_reactants;           // number of reactants
    int         *reactants;             // reactants species
    double      *nu_reactants;          // reactants nus
    double      *ord_reactants;         // reactants orders (nu and FORD/RORD)

    int         *n_products;            // number of products
    int         *products;              // products species
    double      *nu_products;           // products nus
    double      *ord_products;          // products orders (nu and FORD/RORD)

    double      *sum_nu;                // sum of nu from reactants and products

    int         *has_rev_arr;           // flag for reverse arrhenius coefficient
    double      *rev_arr_coeff;         // reverse arrhenius coefficient

    int         *has_high_arr;          // flag for high keyword
    int         *has_low_arr;           // flag for low keyword
    double      *adv_arr_coeff;         // advanced arrhenius coefficients

    int         *has_troe;              // flag for troe coefficients
    int         *n_troe_coeff;          // number of troe coefficients
    double      *troe_coeff;            // Troe coefficients

    int         *has_efficiencies;      // flag for efficiencies
    int         *n_efficiencies;        // number of efficiencies
    int         *sp_efficiencies;       // efficiency species
    double      *efficiencies;          // efficiencies

    int          n_reactions_three;     // Number of three-body reactions
    int         *idx_reactions_three;   // Indices of three-body reactions
    int          n_reactions_low;       // Number of reactions with LOW keyword
    int         *idx_reactions_low;     // Indices of reactions with LOW keyword
    int          n_reactions_high;      // Number of reactions with HIGH keyword
    int         *idx_reactions_high;    // Indices of reactions with HIGH keyword
    int          n_reactions_troe;      // Number of reactions with TROE keyword
    int         *idx_reactions_troe;    // Indices of reactions with TROE keyword
    int          n_reactions_rev;       // Number of reversible reactions
    int         *idx_reactions_rev;     // Indices of reversible reactions
    int          n_reactions_rev_arr;   // Number of reversible reactions with REV keyword
    int         *idx_reactions_rev_arr; // Indices of reversible reactions with REV keyword

    double      *q;                     // rate of progress (mol m-3 s-1)
    double      *k;                     // reaction rates (dep. on reaction)
    double      *pr;                    // reduced pressure
} Reactions_t;

typedef struct Chemistry
{
    Elements_t  *elements;
    Specii_t    *specii;
    Reactions_t *reactions;
} Chemistry_t;

typedef struct State
{
    Chemistry_t *chemistry;     // related chemistry data

    double       p;             // pressure (Pa)
    double       T;             // temperature (K)
    double      *Y;             // mass fractions

    double       R;             // specific gas constant (J kg-1 K-1)
    double       rho;           // density (kg m-3)
    double      *C;             // concentrations (mol m-3)

    double       molar_mass;    // molar mass (kg mol-1)
    double      *X;             // mole fractions ()

    double       cp;            // specific heat capacity (const. pressure, J kg-1 K-1)
    double       cv;            // specific heat capacity (const. volume, J kg-1 K-1)
    double       h;             // specific enthalpy (J kg-1)
    double       s;             // specific entropy (J kg-1 K-1)
    double       u;             // specific inner energy (J kg-1)
    double       g;             // specific free gibbs energy (J kg-1)
} State_t;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
Chemistry_t *read_chemistry_data( const_string_t chemistry_file );
void print_chemistry_info( Chemistry_t *chemistry );
void deallocate_chemistry( Chemistry_t **chemistry );

double calc_sp_cp_r( int i, double T, Chemistry_t *chemistry );
double calc_sp_h_rt( int i, double T, Chemistry_t *chemistry );
double calc_sp_s_r( int i, double T, Chemistry_t *chemistry );
double calc_sp_g_rt( int i, double T, Chemistry_t *chemistry );

void correct_fraction( double *Y, double *Y_out, Chemistry_t *chemistry );
double calc_mix_R( double *Y, Chemistry_t *chemistry );
double calc_mix_molar_mass_from_X( double *X, Chemistry_t *chemistry );
double calc_mix_molar_mass_from_Y( double *Y, Chemistry_t *chemistry );
double calc_mix_cp( double *Y, double T, Chemistry_t *chemistry );
double calc_mix_cv( double cp, double R );
double calc_mix_h( double *Y, double T, Chemistry_t *chemistry );
double calc_mix_s( double *X, double *Y, double p, double T, Chemistry_t *chemistry );
double calc_mix_u( double h, double R, double T );
double calc_mix_g( double h, double s, double T );
void conv_mass_to_mole( double *Y, double mix_molar_mass, double *X, Chemistry_t *chemistry );
void conv_mole_to_mass( double *X, double mix_molar_mass, double *Y, Chemistry_t *chemistry );
void conv_mass_to_conc( double *Y, double rho, double *C, Chemistry_t *chemistry );
void conv_mole_to_conc( double *X, double rho, double molar_mass, double *C, Chemistry_t *chemistry );

void calc_production_rate( double *C, double T, Chemistry_t *chemistry );

double calc_sp_mu( int i, double T, Chemistry_t *chemistry );
double calc_sp_lambda( int i, double p, double T, double mu, double Dii, Chemistry_t *chemistry );
double calc_sp_dii( int i, double p, double T, Chemistry_t *chemistry );
double calc_sp_dij( int i, int j, double p, double T, Chemistry_t *chemistry );

State_t *allocate_state( Chemistry_t *chemistry );
void print_state( State_t *state );
void deallocate_state( State_t **state );
void update_state_isochoric( double p, double rho, double T, double *Y, State_t *state );
void update_state_isobaric( double p, double rho, double T, double *Y, State_t *state );

#endif /* CHEMISTRY_MODULE_H */