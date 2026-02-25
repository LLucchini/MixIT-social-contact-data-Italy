#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "header.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

extern double EPSILON;  
extern int NAGES;
extern gsl_rng * R_GLOBAL;
extern int VERBOSE;
extern int RERUN;
extern double BETA;
extern double REPORTING;
extern double OVERDISP;
extern double SEEDINF;

extern double * PARAMETERS;
extern int N_PARAMETERS;

extern double OMEGA;
extern double GAMMA;
extern double * POP_AGE_INPUT;
extern double * POP_AGE_INPUT_P;
extern double * POP_AGE_INPUT_NP;


extern double POPULATIONSIZE;
extern int NEPICLASSES;
extern double STEP;
extern char * EXP_STRING;

extern int NWEEKS;
//extern int FIT_REPORTING;


extern int WEEK_START_FIT;
extern int WEEK_END_FIT;
extern int NWEEKS_FIT; 
extern int *DATA_INC;
extern unsigned int *CASES_DATA; 

void set_initial_conditions(double *y) {
    int i;
    int offset = NAGES * NEPICLASSES;
    double S_P_total = 0.0;
    double S_NP_total = 0.0;


    // Initialize ODE system
    for(i=0;i<NAGES;i++){ 
      y[0+i*NEPICLASSES]= POP_AGE_INPUT_P[i];                           // S  
      y[1+i*NEPICLASSES]=0.;             // E
      y[2+i*NEPICLASSES]=0.;                           // I
      y[3+i*NEPICLASSES]=0.;                           // R
      y[4+i*NEPICLASSES]=0.;// C 

      y[0+i*NEPICLASSES+NAGES*NEPICLASSES]= POP_AGE_INPUT_NP[i];                           // S
      y[1+i*NEPICLASSES+NAGES*NEPICLASSES]=0.;             // E
      y[2+i*NEPICLASSES+NAGES*NEPICLASSES]=0.;                           // I
      y[3+i*NEPICLASSES+NAGES*NEPICLASSES]=0.;                           // R
      y[4+i*NEPICLASSES+NAGES*NEPICLASSES]=0.;     
    }
     
    for(i = 0; i < NAGES; i++){
        S_P_total  += y[0 + i*NEPICLASSES];
        S_NP_total += y[0 + i*NEPICLASSES + offset];
    }
    
    // Distribute infectious across P and NP
    double N_P = SEEDINF * (S_P_total / (S_P_total + S_NP_total));
    double N_NP = SEEDINF * (S_NP_total / (S_P_total + S_NP_total));

    // Distribute infectious across ages for P
    double S_P_age_total = 0.0;
    for(i = 0; i < NAGES; i++){
      S_P_age_total += y[0 + i*NEPICLASSES];
    }
    
    for(i = 0; i < NAGES; i++){
        double add = 0.0;
        if(S_P_age_total > 0.0){
	  add = N_P * (y[0 + i*NEPICLASSES] / S_P_age_total);
	}
        y[0 + i*NEPICLASSES] -= add;
        y[2 + i*NEPICLASSES] += add; 
    }

    // Distribute infectious across ages for NP
    double S_NP_age_total = 0.0;
    for(i = 0; i < NAGES; i++) S_NP_age_total += y[0 + i*NEPICLASSES + offset];
    for(i = 0; i < NAGES; i++){
        double add = 0.0;
        if(S_NP_age_total > 0.0) add = N_NP * (y[0 + i*NEPICLASSES + offset] / S_NP_age_total);
        y[0 + i*NEPICLASSES + offset] -= add;
        y[2 + i*NEPICLASSES + offset] += add;
    }
    return;

}



void strategy_0(int isim,double *p_model_inc){
  
  char output_file[1000];
  
  register int iage;
  double t, ft, h;
  int nequations,week;
  nequations=NAGES*NEPICLASSES*2;

  /* declarations to solve an ode system*/
  int status;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, nequations);
  gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-6, 0.0);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(nequations);
  
  /* allocate (and initialize) parameters and variables of the system*/
  /*y is the vector of epidemiologica class:  our solution of the ode*/
  double *y;

  double **cumInfectedByAge,  **cumInfectedByAgeP,  **cumInfectedByAgeNP;
  double **susceptibleByAge,**exposedByAge,**infectedByAge,**removedByAge,**populationByAge;
  double *popByAge;
  int countWeekFit=0;
  double cases_model;
  int i;

  y = (double*)calloc(nequations,sizeof(double)); 
  popByAge = (double*)calloc(NAGES,sizeof(double)); 
  
  susceptibleByAge=(double**)calloc(NAGES,sizeof(double*));
  exposedByAge=(double**)calloc(NAGES,sizeof(double*));
  infectedByAge=(double**)calloc(NAGES,sizeof(double*));
  removedByAge=(double**)calloc(NAGES,sizeof(double*));
  cumInfectedByAge=(double**)calloc(NAGES,sizeof(double*));
  cumInfectedByAgeP=(double**)calloc(NAGES,sizeof(double*));
  cumInfectedByAgeNP=(double**)calloc(NAGES,sizeof(double*));
  populationByAge=(double**)calloc(NAGES,sizeof(double*));

  
  for(iage=0;iage<NAGES;iage++){
    susceptibleByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
    exposedByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
    infectedByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
    removedByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
    cumInfectedByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
    cumInfectedByAgeP[iage]=(double*)calloc(NWEEKS,sizeof(double));
    cumInfectedByAgeNP[iage]=(double*)calloc(NWEEKS,sizeof(double));
    populationByAge[iage]=(double*)calloc(NWEEKS,sizeof(double));
  }
  
  for(iage=0;iage<NAGES;iage++){
    for(week=0;week<NWEEKS;week++){
      susceptibleByAge[iage][week]=0.;
      exposedByAge[iage][week]=0.;
      infectedByAge[iage][week]=0.;
      removedByAge[iage][week]=0.;
      cumInfectedByAge[iage][week]=0.;
      cumInfectedByAgeP[iage][week]=0.;
      cumInfectedByAgeNP[iage][week]=0.;
      populationByAge[iage][week]=0.;
    }
  }
 

  /* initialize parameters odes*/
  h=STEP;

  /*Set initial population*/
  set_initial_conditions(y);

  if(VERBOSE==1){

    double tmpP, tmpNP;
    tmpP=0.0;
    tmpNP=0.0;
    
    fprintf(stderr,"check initialization P\n\n");
    
    for(iage=0;iage<NAGES;iage++){
       for(i=0;i<4;i++){
	 fprintf(stderr,"%f\t",y[i+iage*NEPICLASSES]);	 
       }
       tmpP+=y[2+iage*NEPICLASSES];
       fprintf(stderr,"SUM=%f\tPOP=%f",y[0+iage*NEPICLASSES]+y[1+iage*NEPICLASSES]+y[2+iage*NEPICLASSES]+y[3+iage*NEPICLASSES],POP_AGE_INPUT_P[iage]);
       fprintf(stderr,"\n");
    }
   
    fprintf(stderr,"\n\n");
    fprintf(stderr,"check initialization NP\n\n");
    
    for(iage=0;iage<NAGES;iage++){
       for(i=0;i<4;i++){
	 fprintf(stderr,"%f\t",y[i+iage*NEPICLASSES+ NAGES*NEPICLASSES]);
       }
       tmpNP+=y[2+iage*NEPICLASSES+ NAGES*NEPICLASSES];

       fprintf(stderr,"SUM=%f\tPOP=%f",y[0+iage*NEPICLASSES+ NAGES*NEPICLASSES]+y[1+iage*NEPICLASSES+ NAGES*NEPICLASSES]+y[2+iage*NEPICLASSES+ NAGES*NEPICLASSES]+y[3+iage*NEPICLASSES + NAGES*NEPICLASSES],POP_AGE_INPUT_NP[iage]);
       fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n\n");
    fprintf(stderr,"totI0_P=%f\ttotI0_NP=%f\t\n\n",tmpP,tmpNP);
  }

  /*set the system*/
  gsl_odeiv_system sys = {diff_mod, jac, nequations, PARAMETERS};
  t = 0.;

  if(RERUN==0){
    fprintf(stdout,"evaluating isim=%d\n",isim);
    fprintf(stdout,"evaluating BETA=%f\n",BETA); /*lo uso in diff_mod*/
    fprintf(stdout,"evaluating REPORTING=%f\n",REPORTING);
    fprintf(stdout,"evaluating OVERDISP=%f\n",OVERDISP);
    fprintf(stdout,"evaluating SEEDINF=%f\n",SEEDINF);
  }

  /*Print population by age*/
  for(iage=0;iage<NAGES;iage++){
    popByAge[iage]=0.;
  }
  for(iage=0;iage<NAGES;iage++){
    popByAge[iage]=popByAge[iage]+(y[0+iage*NEPICLASSES]+y[1+iage*NEPICLASSES]+y[2+iage*NEPICLASSES]+y[3+iage*NEPICLASSES])+
     (y[0+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[1+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[2+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[3+iage*NEPICLASSES+NAGES*NEPICLASSES]) ;
    if(VERBOSE==1){
      fprintf(stderr,"popByAge[%d]=%f\n", iage,popByAge[iage]);
    }
  }
  
    
  if(VERBOSE==1){
    sprintf(output_file, "%s/output_files/popByAge_CHECK",EXP_STRING);
    fprintf(stderr,"%s\n",output_file);
    save_vector_of_double(output_file,popByAge,NAGES);
  }



  t = 0.0;       
    
  for(week=0;week<NWEEKS;week++){
    for(iage=0;iage<NAGES;iage++){
      y[4+iage*NEPICLASSES]=0.;
      y[4+iage*NEPICLASSES+NAGES*NEPICLASSES]=0.; 
    } 
    t=week;
    ft=week+1;
      while(t < ft){
	status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, ft, &h, y);
	if (status != GSL_SUCCESS) break;
	h = STEP; 
      }

    for(iage=0;iage<NAGES;iage++){
    
      cumInfectedByAgeP[iage][week]=y[4+iage*NEPICLASSES];
      cumInfectedByAgeNP[iage][week]=y[4+iage*NEPICLASSES+NAGES*NEPICLASSES];
      cumInfectedByAge[iage][week]= cumInfectedByAgeP[iage][week]+ cumInfectedByAgeNP[iage][week];
    }
    
    
    for(iage=0;iage<NAGES;iage++){
      susceptibleByAge[iage][week]=y[0+iage*NEPICLASSES]+y[0+iage*NEPICLASSES+NAGES*NEPICLASSES];
      exposedByAge[iage][week]=y[1+iage*NEPICLASSES]+y[1+iage*NEPICLASSES+NAGES*NEPICLASSES];
      infectedByAge[iage][week]=y[2+iage*NEPICLASSES]+y[2+iage*NEPICLASSES+NAGES*NEPICLASSES];
      removedByAge[iage][week]=y[3+iage*NEPICLASSES]+y[3+iage*NEPICLASSES+NAGES*NEPICLASSES];
      populationByAge[iage][week]=(y[0+iage*NEPICLASSES]+y[1+iage*NEPICLASSES]+y[2+iage*NEPICLASSES]+y[3+iage*NEPICLASSES])+(y[0+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[1+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[2+iage*NEPICLASSES+NAGES*NEPICLASSES]+y[3+iage*NEPICLASSES+NAGES*NEPICLASSES]);
      if (fabs(populationByAge[iage][week] - POP_AGE_INPUT[iage]) > EPSILON) {
	fprintf(stderr, "populationByAge[%d][%d]=%.12f\t", iage, week, populationByAge[iage][week]);
	fprintf(stderr, "POPINPUT=%.12f\n", POP_AGE_INPUT[iage]);
	fprintf(stderr, "ERRORE!!\n");
	exit(1);
      }
    }

    if(week>=WEEK_START_FIT && week<=WEEK_END_FIT){
	cases_model=0.;
	for(i=0;i<NAGES;i++){
	  if(cumInfectedByAge[i][week]>0)
	    cases_model+=REPORTING*cumInfectedByAge[i][week];
	}
	if(VERBOSE==1){
	  fprintf(stderr,"week=%d\tcases_data=%u\tall_cases_model=%f\treported_cases_model=%f\n",week,CASES_DATA[week],cases_model/REPORTING,cases_model);
	}
	
	

	if(RERUN==0){
	  p_model_inc[countWeekFit]=OVERDISP/(cases_model+OVERDISP);
	}
	if(VERBOSE==1){
	  fprintf(stderr,"p_model_inc[%d]=%f\n",countWeekFit,p_model_inc[countWeekFit]);
	}

      
	countWeekFit++;  
    }
  
 
  }

  /*PRINT OUTPUT_FILES*/
  sprintf(output_file, "%s/output_files/CumInfByAgeP_%d",EXP_STRING,isim);
  save_matrix_of_double(output_file,cumInfectedByAgeP,NAGES,NWEEKS);

  sprintf(output_file, "%s/output_files/CumInfByAgeNP_%d",EXP_STRING,isim);
  save_matrix_of_double(output_file,cumInfectedByAgeNP,NAGES,NWEEKS);
     
  if(VERBOSE==1){
    sprintf(output_file, "%s/output_files/susceptibleByAge_%d",EXP_STRING,isim);
    save_matrix_of_double(output_file,susceptibleByAge,NAGES,NWEEKS);

    sprintf(output_file, "%s/output_files/exposeddByAge_%d",EXP_STRING,isim);
    save_matrix_of_double(output_file,exposedByAge,NAGES,NWEEKS);

    sprintf(output_file, "%s/output_files/infectedByAge_%d",EXP_STRING,isim);
    save_matrix_of_double(output_file,infectedByAge,NAGES,NWEEKS);
  
    sprintf(output_file, "%s/output_files/removedByAge_%d",EXP_STRING,isim);
    save_matrix_of_double(output_file,removedByAge,NAGES,NWEEKS);
  }
  
  if(isim==0){
    sprintf(output_file, "%s/output_files/populationByAge_%d",EXP_STRING,isim);
    save_matrix_of_double(output_file,populationByAge,NAGES,NWEEKS);
  }
  


  fflush(stderr);
  fflush(stdout);


 
  free(y);
  free(popByAge);

  
  for(iage=0;iage<NAGES;iage++){
    free(cumInfectedByAge[iage]);
    free(susceptibleByAge[iage]);
    free(exposedByAge[iage]);
    free(infectedByAge[iage]);
    free(removedByAge[iage]);
    free(populationByAge[iage]);
  }
  
  free(cumInfectedByAge);
  free(susceptibleByAge);
  free(exposedByAge);
  free(infectedByAge);
  free(removedByAge);
  free(populationByAge);


  

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);  
  
  return;
}





int jac(double t, const double y[], double *dfdy, double dfdt[], void *p){
  (void)t;
  (void)y;
  (void)p;

  int nequations;
  nequations=NAGES*NEPICLASSES;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, nequations, nequations);
  gsl_matrix * m = &dfdy_mat.matrix;
  register unsigned short int i,j;
  for(i=0; i<nequations; i++)
    for(j=0; j<nequations; j++)
      gsl_matrix_set(m, i, j, 0.);
  for(i=0; i<nequations; i++)
    dfdt[i] = 0.;
  return GSL_SUCCESS;
} 
 
