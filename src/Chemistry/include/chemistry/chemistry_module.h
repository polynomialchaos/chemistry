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
#define RM NA *KB         /** universal gas constant, mol-1 m2 kg s-2 K-1 */
#define P0 101325.0       /** standard pressure, Pa */
#define T0 273.15         /** standard temperature, K */

#define BOUNDS 3 /** number of temperature bounds */
#define NASA 7   /** number of NASA polynomial coefficients */
#define ARR 3    /** number of Arrhenius coefficients */
#define TROE 5   /** number of Troe pressure function coefficients */
#define KINDIM 2 /** nmuber of kinetic dimensions (forward/backward) */

#define YONE 1e-6     /** maximum tolerance for sum y - 1 */
#define YSMALL 1e-100 /** minimum species concentration */

/*******************************************************************************
 * @brief Elements structure
 ******************************************************************************/
typedef struct Elements
{
    int n_elements; /** number of elements */

    string_t *symbol; /** symbol */
    double *mass;     /** mass, kg mol-1 */
} elements_t;

/*******************************************************************************
 * @brief Reactions structure
 ******************************************************************************/
typedef struct Reactions
{
    int n_reactions;      /** number of reactions */
    int max_reactants;    /** maximum number of reactants */
    int max_products;     /** maximum number of products */
    int max_troe_coeff;   /** maximum number of troe coefficients */
    int max_efficiencies; /** maximum number of efficiencies */

    int *type;            /** reaction type */
    int *falloff_species; /** fall off species */

    double *arr_coeff;  /** Arrhenius coefficient */
    int *is_reversible; /** flag for reversibility */

    int *n_reactants;      /** number of reactants */
    int *reactants;        /** reactants species */
    double *nu_reactants;  /** reactants nus */
    double *ord_reactants; /** reactants orders (nu and FORD/RORD) */

    int *n_products;      /** number of products */
    int *products;        /** products species */
    double *nu_products;  /** products nus */
    double *ord_products; /** products orders (nu and FORD/RORD) */

    double *sum_nu; /** sum of nu from reactants and products */

    int *has_rev_arr;      /** flag for reverse arrhenius coefficient */
    double *rev_arr_coeff; /** reverse arrhenius coefficient */

    int *has_high_arr;     /** flag for high keyword */
    int *has_low_arr;      /** flag for low keyword */
    double *adv_arr_coeff; /** advanced arrhenius coefficients */

    int *has_troe;      /** flag for troe coefficients */
    int *n_troe_coeff;  /** number of troe coefficients */
    double *troe_coeff; /** Troe coefficients */

    int *has_efficiencies; /** flag for efficiencies */
    int *n_efficiencies;   /** number of efficiencies */
    int *sp_efficiencies;  /** efficiency species */
    double *efficiencies;  /** efficiencies */

    int n_reactions_three;      /** Number of three-body reactions */
    int *idx_reactions_three;   /** Indices of three-body reactions */
    int n_reactions_low;        /** Number of reactions with LOW keyword */
    int *idx_reactions_low;     /** Indices of reactions with LOW keyword */
    int n_reactions_high;       /** Number of reactions with HIGH keyword */
    int *idx_reactions_high;    /** Indices of reactions with HIGH keyword */
    int n_reactions_troe;       /** Number of reactions with TROE keyword */
    int *idx_reactions_troe;    /** Indices of reactions with TROE keyword */
    int n_reactions_rev;        /** Number of reversible reactions */
    int *idx_reactions_rev;     /** Indices of reversible reactions */
    int n_reactions_rev_arr;    /** Number of reversible reactions with REV keyword */
    int *idx_reactions_rev_arr; /** Indices of reversible reactions with REV keyword */

    double *q;  /** rate of progress (mol m-3 s-1) */
    double *k;  /** reaction rates (dep. on reaction) */
    double *pr; /** reduced pressure */
} reactions_t;

/*******************************************************************************
 * @brief Specii structure
 ******************************************************************************/
typedef struct Specii
{
    int n_specii;        /** number of specii */
    int max_reac_points; /** maximum number of reaction points */

    string_t *symbol; /** symbol */
    int *is_inert;    /** flag for inert species */

    int *composition;   /** elementary composition */
    int *phase;         /** phase */
    double *bounds;     /** temperature bounds for NASA polynomials */
    double *coeff_high; /** coefficients for high temperature range */
    double *coeff_low;  /** coefficients for low temperature range */

    int *geom;       /** geometrical configuration */
    double *pot_lj;  /** lennard-Jones potential */
    double *col_lj;  /** lennard-Jones collision diameter */
    double *dip_mo;  /** dipole momentum */
    double *pol;     /** polarizeability */
    double *rot_rel; /** rotational relaxation */

    int *n_reac_points;     /** number of reaction points */
    int *reac_points;       /** reaction points */
    double *nu_reac_points; /** reaction points nu */

    double *molar_mass;      /** molar mass of species, kg mol-1 */
    double *molecule_weight; /** molecular mass of species, kg */
    double *Rsp;             /** specific gas constant, J kg-1 K-1 */

    double *omega; /** production rates (kg m-3 s-1) */
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
    chemistry_t *chemistry; /** related chemistry data */

    double p;  /** pressure (Pa) */
    double T;  /** temperature (K) */
    double *Y; /** mass fractions */

    double R;   /** specific gas constant (J kg-1 K-1) */
    double rho; /** density (kg m-3) */
    double *C;  /** concentrations (mol m-3) */

    double molar_mass; /** molar mass (kg mol-1) */
    double *X;         /** mole fractions () */

    double cp; /** specific heat capacity (const. pressure, J kg-1 K-1) */
    double cv; /** specific heat capacity (const. volume, J kg-1 K-1) */
    double h;  /** specific enthalpy (J kg-1) */
    double s;  /** specific entropy (J kg-1 K-1) */
    double u;  /** specific inner energy (J kg-1) */
    double g;  /** specific free gibbs energy (J kg-1) */
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
 * @brief Return dimensionless enthalpy
 * @param i
 * @param T
 * @param chemistry
 * @return double
 ******************************************************************************/
double calc_species_h_rt(int i, double T, chemistry_t *chemistry);

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
double calc_sp_g_rt(int i, double T, chemistry_t *chemistry);

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

double calc_sp_mu(int i, double T, chemistry_t *chemistry);
double calc_sp_lambda(int i, double p, double T, double mu, double Dii, chemistry_t *chemistry);
double calc_sp_dii(int i, double p, double T, chemistry_t *chemistry);
double calc_sp_dij(int i, int j, double p, double T, chemistry_t *chemistry);

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