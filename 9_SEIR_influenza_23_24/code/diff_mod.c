#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv.h>
#include "header.h"

extern double POPULATIONSIZE;
extern double * POP_AGE_INPUT;
extern double * POP_AGE_INPUT_P;
extern double * POP_AGE_INPUT_NP;
extern double * SUSC_AGE_INPUT;


extern double EPSILON;
extern double *  PARAMETERS;
extern int NAGES;
extern int NEPICLASSES;
extern double BETA;
extern double OMEGA;
extern double GAMMA;
extern double **C_MATRIX_PRESENCE;
extern double **C_MATRIX_NOT_PRESENCE;


int diff_mod(double t, const double y[], double f[], void *params){
    (void) t;
    (void) params;

    double *totp_P = (double*)calloc(NAGES, sizeof(double));
    double *totp_NP = (double*)calloc(NAGES, sizeof(double));
    double *totp = (double*)calloc(NAGES, sizeof(double));
    
    double *foi_P  = (double*)calloc(NAGES, sizeof(double));
    double *foi_NP = (double*)calloc(NAGES, sizeof(double));

    int i, j;
    double tmp;

    for(i = 0; i < NAGES; i++){
      totp_P[i]  = y[0 + i*NEPICLASSES] + y[1 + i*NEPICLASSES] + y[2 + i*NEPICLASSES] + y[3 + i*NEPICLASSES];
      totp_NP[i] = y[0 + i*NEPICLASSES + NAGES*NEPICLASSES] + 
	y[1 + i*NEPICLASSES + NAGES*NEPICLASSES] +
	y[2 + i*NEPICLASSES + NAGES*NEPICLASSES] +
	y[3 + i*NEPICLASSES + NAGES*NEPICLASSES];
      totp[i]=totp_P[i]+totp_NP[i];
      if (fabs(totp_P[i] - POP_AGE_INPUT_P[i]) > EPSILON) {
	fprintf(stderr, "totp_P[%d]=%.12f\t", i, totp_P[i]);
	fprintf(stderr, "POPINPUT_P=%.12f\n", POP_AGE_INPUT_P[i]);
	fprintf(stderr, "ERROR POP DIFF MOD!!\n");
	exit(1);
      }
      if (fabs(totp_NP[i] - POP_AGE_INPUT_NP[i]) > EPSILON) {
	fprintf(stderr, "totp_NP[%d]=%.12f\t", i, totp_NP[i]);
	fprintf(stderr, "POPINPUT_NP=%.12f\n", POP_AGE_INPUT_NP[i]);
	fprintf(stderr, "ERROR POP DIFF MOD!!\n");
	exit(1);
      }
      if (fabs(totp[i] - POP_AGE_INPUT[i]) > EPSILON) {
	fprintf(stderr, "totp[%d]=%.12f\t", i, totp[i]);
	fprintf(stderr, "POPINPUT=%.12f\n", POP_AGE_INPUT[i]);
	fprintf(stderr, "ERROR POP DIFF MOD!!\n");
	exit(1);
      }

    }


    // Compute FOI for group P
    for(i = 0; i < NAGES; i++){
        tmp = 0.0;
        for(j = 0; j < NAGES; j++){
	  if(totp[j] > 0.0){
	    tmp += SUSC_AGE_INPUT[i]*C_MATRIX_PRESENCE[i][j] * ((y[2 + j*NEPICLASSES]+ y[2 + j*NEPICLASSES + NAGES*NEPICLASSES])/totp[j]);
	  }
        }
        foi_P[i] = BETA * tmp;
    }

     // Compute FOI for group NP
    for(i = 0; i < NAGES; i++){
        tmp = 0.0;
        for(j = 0; j < NAGES; j++){
	  if(totp[j] > 0.0){
	    tmp += SUSC_AGE_INPUT[i]* C_MATRIX_NOT_PRESENCE[i][j] * ((y[2 + j*NEPICLASSES]+ y[2 + j*NEPICLASSES + NAGES*NEPICLASSES])/totp[j]);
	  }
	}
        foi_NP[i] = BETA * tmp;
    }

    // ODE system for group P
    for(i = 0; i < NAGES; i++){
        f[0 + i*NEPICLASSES] = -foi_P[i] * y[0 + i*NEPICLASSES];                      // S
        f[1 + i*NEPICLASSES] =  foi_P[i] * y[0 + i*NEPICLASSES] - OMEGA * y[1 + i*NEPICLASSES]; // E
        f[2 + i*NEPICLASSES] =  OMEGA * y[1 + i*NEPICLASSES] - GAMMA * y[2 + i*NEPICLASSES];  // I
        f[3 + i*NEPICLASSES] =  GAMMA * y[2 + i*NEPICLASSES];                        // R
        f[4 + i*NEPICLASSES] =  foi_P[i] * y[0 + i*NEPICLASSES];                     // incidence
    }

    // ODE system for group NP
    for(i = 0; i < NAGES; i++){
        int offset = NAGES*NEPICLASSES;
        f[0 + i*NEPICLASSES + offset] = -foi_NP[i] * y[0 + i*NEPICLASSES + offset];                      // S
        f[1 + i*NEPICLASSES + offset] =  foi_NP[i] * y[0 + i*NEPICLASSES + offset] - OMEGA * y[1 + i*NEPICLASSES + offset]; // E
        f[2 + i*NEPICLASSES + offset] =  OMEGA * y[1 + i*NEPICLASSES + offset] - GAMMA * y[2 + i*NEPICLASSES + offset];  // I
        f[3 + i*NEPICLASSES + offset] =  GAMMA * y[2 + i*NEPICLASSES + offset];                        // R
        f[4 + i*NEPICLASSES + offset] =  foi_NP[i] * y[0 + i*NEPICLASSES + offset];                     // incidence
    }

    free(totp_P); free(totp_NP); free(totp);
    free(foi_P); free(foi_NP);

    return GSL_SUCCESS;
}
