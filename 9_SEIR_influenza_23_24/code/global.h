/*General variables*/
gsl_rng * R_GLOBAL;               // Global random number generator (from the GNU Scientific Library)
int ISEED;                        // Seed for the random number generator
int IEXP;                         // Index or ID of the experiment
char * EXP_STRING;                // String label describing the experiment
int VERBOSE;                      // =1 to print diagnostic checks, =0 to suppress them
//int FIT_REPORTING;                // Controls whether to include reporting in model fitting
int RERUN;                        // Flag to ontrols whether to perform model calibration (0) or run simulations using known parameters (!=0).

/*Contact matrices*/
double ** C_MATRIX_PRESENCE;      // Contact matrix for in-person population group
double ** C_MATRIX_NOT_PRESENCE;  // Contact matrix for not in-person population group

/*Simulation parameters*/
int NSIM;                         // Number of MCMC (Markov Chain Monte Carlo) simulations

/*Time variables*/
int NWEEKS = 24;                  // Total number of weeks simulated
int NWEEKS_FIT = 24;              // Number of weeks used for model fitting (weeks 46–52 of 2023 and 1–17 of 2024)
int WEEK_START_FIT = 0;           // Starting week index for fitting
int WEEK_END_FIT = 23;            // Ending week index for fitting

/*Epidemiological variables*/
double OMEGA;                     // Latency rate (in weeks^-1)
double GAMMA;                     // Recovery rate (in weeks^-1)
int NEPICLASSES = 5;              // Number of epidemiological compartments (SEIRC model)
                                  // C = cumulative weekly infections

/*Population variables*/
double POPULATIONSIZE;            // Total population size
double* POP_AGE_INPUT;            // Age-structured population input
double* POP_AGE_INPUT_P;          // Population in-person by age group
double* POP_AGE_INPUT_NP;         // Population not in-person by age group
//double* REPORTING_AGE_INPUT;      // Age-specific reporting probabilities
double* SUSC_AGE_INPUT;           // Age-specific susceptibility values

int NAGES = 15;                   // Number of age groups (0–4y to 65–69y, plus 70+)

/*ODE system*/
double STEP = 0.1;                // Time step for ODE solver (in weeks; e.g., 1/7 for daily)
double EPSILON = 1e-2;            // Tolerance for population conservation errors

/*MCMC variables*/
int IDCURRENT;                    // Current MCMC iteration index
long double LKHCURRENT;           // Current log-likelihood value
long double LKHCANDIDATE;         // Candidate log-likelihood value
double BETA;                      // Scalinf factor for transmissibility
double BETACURRENT;               // Current value of BETA
double BETACANDIDATE;             // Proposed candidate for BETA
double REPORTING;                 // Reporting ratio parameter
double REPORTINGCURRENT;          // Current value of REPORTING
double REPORTINGCANDIDATE;        // Proposed candidate for REPORTING
double OVERDISP;                  // Overdispersion parameter (for likelihood)
double OVERDISPCURRENT;           // Current value of OVERDISPERSION
double OVERDISPCANDIDATE;         // Proposed candidate for OVERDISPERSION
double SEEDINF;                   // Infectious individuals at simulation start
double SEEDINFCURRENT;            // Current value of SEEDINF
double SEEDINFCANDIDATE;          // Proposed candidate for SEEDINF
double BETASTART;                 // Initial starting value for BETA
double REPORTINGSTART;            // Initial starting value for REPORTING
double OVERDISPSTART;             // Initial starting value for OVERDISP
double SEEDINFSTART;              // Initial starting value for SEEDINF
double SD_BETA;                   // Proposal standard deviation for BETA
double SD_REPORTING;              // Proposal standard deviation for REPORTING
double SD_OVERDISP;               // Proposal standard deviation for OVERDISP
double SD_SEEDINF;                // Proposal standard deviation for SEEDINF

double * PARAMETERS;              // Vector storing model parameters
int NPARAMETERS = 4;              // Number of model parameters

/*Epidemiological data input*/
double *DATA_INC;                 // Observed weekly incidence data (ILI+_A for 2023–2024)
unsigned int *CASES_DATA;         // case count data
