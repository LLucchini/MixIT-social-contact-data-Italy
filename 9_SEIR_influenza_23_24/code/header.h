#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

int parser(int n,char *array[],char ***flag,char ***value,int *nflags);
char *get_value(char *flag[],char *value[],int nflags,char opt[]);
void percent(int n,int d,int s, FILE *out);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *p);
int diff_mod(double t, const double y[], double f[], void *params);
void mcmc();
void rerun();
void save_vector_of_double(const char *filename,const double *vector,unsigned maxlen);
void save_matrix_of_double(const char *filename, double **matrix,unsigned nrow, unsigned ncol);
void save_vector_of_int(const char *filename,const int *vector,unsigned maxlen);
void read_vector_of_double(double * vector, int dimvector, char * filename); 
void read_vector_of_int(int * vector, int dimvector, char * filename);
void strategy_0(int isim, double *p_model_inc);
int read_demo_data(char sep);
int read_epi_data(char sep);
int read_incidence();

/* int read_sero_data_varicella(char sep); */
int get_line(char **line,FILE *fp);
void update_c_matrix();
int read_param(char file[],char sep);
int readMatrix(const char *filename,
               double *beta,
               double *reporting,
               double *seedinf,
               size_t nrows);
