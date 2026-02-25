#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <math.h>
#include "header.h"
extern double POPULATIONSIZE;
extern int NSIM;

extern double BETA;
extern double REPORTING;
extern double SEEDINF;

extern gsl_rng * R_GLOBAL;
extern int ISEED;
extern int IEXP;
extern char * EXP_STRING;
extern double OMEGA;
extern double GAMMA;

extern int VERBOSE;
extern int NWEEKS_FIT;




void rerun(){


  char * directory_output;
  directory_output = (char*)calloc(1000,sizeof(char)); 

  double *beta;
  double *reporting;
  double *seedinf;
  double *p_model_inc;
  int isim;
 
  beta=(double *)calloc(NSIM, sizeof(double));
  reporting=(double *)calloc(NSIM, sizeof(double));
  seedinf=(double *)calloc(NSIM, sizeof(double));
  p_model_inc=(double *)calloc(NWEEKS_FIT, sizeof(double));

  sprintf(directory_output, "%soutput_files",EXP_STRING);
  fprintf(stderr,"%s\n",directory_output);
  if(mkdir(directory_output,S_IRWXU) != 0 && errno != EEXIST){
    perror("mkdir failed");
    exit(1);
  }
  if (readMatrix("input_files/mcmcrerun", beta, reporting, seedinf, NSIM) == 0) {
    printf("File read successfully.\n");
    printf("First beta value: %f\n", beta[0]);
    printf("First reporting value: %f\n", reporting[0]);
    printf("First seedinf value: %f\n", seedinf[0]);
  }
  for(isim=0;isim<NSIM;isim++){
    if(VERBOSE==1){
      printf("isim=%d\tbeta=%f\t;reporting=%f\tseedinf=%f\n",isim,beta[isim], reporting[isim], seedinf[isim]);
    }
    BETA=beta[isim];
    REPORTING=reporting[isim];
    SEEDINF=seedinf[isim];
    strategy_0(isim,p_model_inc);

  }  
  // Free memory
  free(beta);
  free(reporting);
  free(seedinf);

  free(p_model_inc);
    
  
  free(directory_output);
  return;
}
