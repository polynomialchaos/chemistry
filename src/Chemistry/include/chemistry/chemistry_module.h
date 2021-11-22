/*******************************************************************************
 * @file chemistry_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef CHEMISTRY_MODULE_H
#define CHEMISTRY_MODULE_H

#include "basec/basec_module.h"

#define NA 6.022140857e23 /** Avogadro number, mol-1 */
#define KB 1.38064852e-23 /** Boltzmann constant, m2 kg s-2 K-1 */
#define RM NA *KB         /** Universal gas constant, mol-1 m2 kg s-2 K-1 */
#define P0 101325.0       /** Standard pressure, Pa */
#define T0 273.15         /** Standard temperature, K */

#define BOUNDS 3 /** Number of temperature bounds */
#define NASA 7   /** Number of NASA polynomial coefficients */
#define ARR 3    /** Number of Arrhenius coefficients */
#define TROE 5   /** Number of Troe pressure function coefficients */
#define KINDIM 2 /** Nmuber of kinetic dimensions (forward/backward) */

#define YONE 1e-6     /** Maximum tolerance for sum y - 1 */
#define YSMALL 1e-100 /** Minimum species concentration */

/*******************************************************************************
 * @brief Elements structure
 ******************************************************************************/
typedef struct Elements
{
    int n_elements;   /** Number of elements */
    string_t *symbol; /** Symbol */
    double *mass;     /** Mass, kg mol-1 */
} elements_t;

/*******************************************************************************
 * @brief Reactions structure
 ******************************************************************************/
typedef struct Reactions
{
    int n_reactions;          /** Number of reactions */
    int max_reactants;        /** Maximum number of reactants */
    int max_products;         /** Maximum number of products */
    int max_troe_coeff;       /** Maximum number of troe coefficients */
    int max_efficiencies;     /** Maximum number of efficiencies */
    int *type;                /** Reaction type */
    int *falloff_species;     /** Fall off species */
    double *arr_coeff;        /** Arrhenius coefficient */
    int *is_reversible;       /** Flag for reversibility */
    int *n_reactants;         /** Number of reactants */
    int *reactants;           /** Reactants species */
    double *nu_reactants;     /** Reactants nus */
    double *ord_reactants;    /** Reactants orders (nu and FORD/RORD) */
    int *n_products;          /** Number of products */
    int *products;            /** Products species */
    double *nu_products;      /** Products nus */
    double *ord_products;     /** Products orders (nu and FORD/RORD) */
    double *sum_nu;           /** Sum of nu from reactants and products */
    int *has_rev_arr;         /** Flag for reverse arrhenius coefficient */
    double *rev_arr_coeff;    /** Reverse arrhenius coefficient */
    int *has_high_arr;        /** Flag for high keyword */
    int *has_low_arr;         /** Flag for low keyword */
    double *adv_arr_coeff;    /** Advanced arrhenius coefficients */
    int *has_troe;            /** Flag for troe coefficients */
    int *n_troe_coeff;        /** Number of troe coefficients */
    double *troe_coeff;       /** Troe coefficients */
    int *has_efficiencies;    /** Flag for efficiencies */
    int *n_efficiencies;      /** Number of efficiencies */
    int *sp_efficiencies;     /** Efficiency species */
    double *efficiencies;     /** Efficiencies */
    int n_reactions_three;    /** Number of three-body reactions */
    int *idx_reactions_three; /** Three-body reactions */
    int n_reactions_low;      /** Number of reactions with LOW keyword */
    int *idx_reactions_low;   /** Reactions with LOW keyword */
    int n_reactions_high;     /** Number of reactions with HIGH keyword */
    int *idx_reactions_high;  /** Reactions with HIGH keyword */
    int n_reactions_troe;     /** Number of reactions with TROE keyword */
    int *idx_reactions_troe;  /** Reactions with TROE keyword */
    int n_reactions_rev;      /** Number of reversible reactions */
    int *idx_reactions_rev;   /** Reversible reactions */
    /** Number of reversible reactions with REV keyword */
    int n_reactions_rev_arr;
    int *idx_reactions_rev_arr; /** Reversible reactions with REV keyword */
    double *q;                  /** Rate of progress, mol m-3 s-1 */
    double *k;                  /** Reaction rates, dep. on reaction */
    double *pr;                 /** Reduced pressure */
} reactions_t;

/*******************************************************************************
 * @brief Specii structure
 ******************************************************************************/
typedef struct Specii
{
    int n_specii;            /** Number of specii */
    int max_reac_points;     /** Maximum number of reaction points */
    string_t *symbol;        /** Symbol */
    int *is_inert;           /** Flag for inert species */
    int *composition;        /** Elementary composition */
    int *phase;              /** Phase */
    double *bounds;          /** Temperature bounds for NASA polynomials */
    double *coeff_high;      /** Coefficients for high temperature range */
    double *coeff_low;       /** Coefficients for low temperature range */
    int *geom;               /** Geometrical configuration */
    double *pot_lj;          /** Lennard-Jones potential */
    double *col_lj;          /** Lennard-Jones collision diameter */
    double *dip_mo;          /** Dipole momentum */
    double *pol;             /** Polarizeability */
    double *rot_rel;         /** Rotational relaxation */
    int *n_reac_points;      /** Number of reaction points */
    int *reac_points;        /** Reaction points */
    double *nu_reac_points;  /** Reaction points nu */
    double *molar_mass;      /** Molar mass of species, kg mol-1 */
    double *molecule_weight; /** Molecular mass of species, kg */
    double *Rsp;             /** Specific gas constant, J kg-1 K-1 */
    double *omega;           /** Production rates, kg m-3 s-1 */
} specii_t;

/*******************************************************************************
 * @brief Chemistry structure
 ******************************************************************************/
typedef struct Chemistry
{
    elements_t *elements;
    specii_t *specii;
    reactions_t *reactions;
} chemistry_t;

/*******************************************************************************
 * @brief State structure
 ******************************************************************************/
typedef struct State
{
    chemistry_t *chemistry; /** Related chemistry data */
    double p;               /** Pressure, Pa */
    double T;               /** Temperature, K */
    double *Y;              /** Mass fractions */
    double R;               /** Specific gas constant, J kg-1 K-1 */
    double rho;             /** Density, kg m-3 */
    double *C;              /** Concentrations, mol m-3 */
    double molar_mass;      /** Molar mass, kg mol-1 */
    double *X;              /** Mole fractions */
    /** Specific heat capacity (const. pressure), J kg-1 K-1 */
    double cp;
    /** Specific heat capacity (const. volume), J kg-1 K-1 */
    double cv;
    double h; /** Specific enthalpy, J kg-1 */
    double s; /** Specific entropy, J kg-1 K-1 */
    double u; /** Specific inner energy, J kg-1 */
    double g; /** Specific free gibbs energy, J kg-1 */
} state_t;

/*******************************************************************************
 * @brief Allocate the state structure
 * @param chemistry
 * @return state_t*
 ******************************************************************************/
state_t *allocate_state(chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures heat capacity (const. pressure)
 * @param Y
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_cp(double *Y, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures heat capacity (const. volume)
 * @param cp
 * @param R
 * @return double
 ******************************************************************************/
double calc_mix_cv(double cp, double R);

/*******************************************************************************
 * @brief Calculate the mixtures free Gibb's energy
 * @param h
 * @param s
 * @param T
 * @return double
 ******************************************************************************/
double calc_mix_g(double h, double s, double T);

/*******************************************************************************
 * @brief Calculate the mixtures enthalpy
 * @param Y
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_h(double *Y, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures molar mass from mole fraction
 * @param X
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_molar_mass_from_X(double *X, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures molar mass from mass fraction
 * @param Y
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_molar_mass_from_Y(double *Y, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures gas constant
 * @param Y
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_R(double *Y, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures entropy (with Dalton's fix)
 * @param X
 * @param Y
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_mix_s(double *X, double *Y,
                  double p, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate the mixtures internal energy
 * @param h
 * @param R
 * @param T
 * @return double
 ******************************************************************************/
double calc_mix_u(double h, double R, double T);

/*******************************************************************************
 * @brief Calculate the production rates
 * @param C
 * @param T
 * @param chemistry
 ******************************************************************************/
void calc_production_rate(double *C, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Return dimensionless heat capacity
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_cp_r(int i, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate species self diffusion coefficient
 * @param i
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_dii(int i, double p, double T,
                        chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate species binary diffusion coefficient
 * @param i
 * @param j
 * @param p
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_dij(int i, int j, double p, double T,
                        chemistry_t *chemistry);

/*******************************************************************************
 * @brief Return dimensionless enthalpy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_h_rt(int i, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate species heat conductivity
 * @param i
 * @param p
 * @param T
 * @param mu
 * @param Dii
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_lambda(int i, double p, double T, double mu, double Dii,
                           chemistry_t *chemistry);

/*******************************************************************************
 * @brief Calculate species viscosity
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_mu(int i, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Return dimensionless entropy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_s_r(int i, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Return dimensionless free gibb's energy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_g_rt(int i, double T, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Convert the give mass fractions to concentrations
 * @param Y
 * @param rho
 * @param C
 * @param chemistry
 ******************************************************************************/
void conv_mass_to_conc(double *Y, double rho,
                       double *C, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Convert the give mass fractions to mole fractions
 * @param Y
 * @param mix_molar_mass
 * @param X
 * @param chemistry
 ******************************************************************************/
void conv_mass_to_mole(double *Y, double mix_molar_mass,
                       double *X, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Convert the give mole fractions to concentrations
 * @param X
 * @param rho
 * @param molar_mass
 * @param C
 * @param chemistry
 ******************************************************************************/
void conv_mole_to_conc(double *X, double rho, double molar_mass,
                       double *C, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Convert the give mole fractions to concentrations
 * @param X
 * @param mix_molar_mass
 * @param Y
 * @param chemistry
 ******************************************************************************/
void conv_mole_to_mass(double *X, double mix_molar_mass,
                       double *Y, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Correct the given mass/mole fractions
 * @param Y
 * @param Y_out
 * @param chemistry
 ******************************************************************************/
void correct_fraction(double *Y, double *Y_out, chemistry_t *chemistry);

/*******************************************************************************
 * @brief Deallocate chemistry structure
 * @param chemistry
 ******************************************************************************/
void deallocate_chemistry(chemistry_t *chemistry);

/*******************************************************************************
 * @brief Deallocate the state structure
 * @param state
 ******************************************************************************/
void deallocate_state(state_t *state);

/*******************************************************************************
 * @brief Print chemistry structure
 * @param chemistry
 ******************************************************************************/
void print_chemistry_info(chemistry_t *chemistry);

/*******************************************************************************
 * @brief Print the state structure
 * @param state
 ******************************************************************************/
void print_state(state_t *state);

/*******************************************************************************
 * @brief Read chemistry file and fill chemistry structure
 * @param chemistry_file
 * @return chemistry_t*
 ******************************************************************************/
chemistry_t *read_chemistry_data(cstring_t chemistry_file);

/*******************************************************************************
 * @brief Update the state structure (isochoric)
 * @param p
 * @param rho
 * @param T
 * @param Y
 * @param state
 ******************************************************************************/
void update_state_isochoric(double p, double rho, double T,
                            double *Y, state_t *state);

/*******************************************************************************
 * @brief Update the state structure (isobaric)
 * @param p
 * @param rho
 * @param T
 * @param Y
 * @param state
 ******************************************************************************/
void update_state_isobaric(double p, double rho, double T,
                           double *Y, state_t *state);

#endif /* CHEMISTRY_MODULE_H */