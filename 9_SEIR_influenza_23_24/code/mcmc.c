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
extern double OVERDISP;
extern double SEEDINF;
extern double BETACURRENT;
extern double REPORTINGCURRENT;
extern double OVERDISPCURRENT;
extern double SEEDINFCURRENT;
extern double BETACANDIDATE;
extern double REPORTINGCANDIDATE;
extern double OVERDISPCANDIDATE;
extern double SEEDINFCANDIDATE;
extern double BETASTART;
extern double REPORTINGSTART;
extern double OVERDISPSTART;
extern double SEEDINFSTART;

extern double SD_BETA;
extern double SD_REPORTING;
extern double SD_OVERDISP;
extern double SD_SEEDINF;

extern int IDCURRENT;
extern long double LKHCURRENT;
extern long double LKHCANDIDATE;
extern gsl_rng * R_GLOBAL;
extern int ISEED;
extern int IEXP;
extern char * EXP_STRING;
extern double OMEGA;
extern double GAMMA;

extern int VERBOSE;
extern int NWEEKS_FIT;
extern double *DATA_INC; 
extern unsigned int *CASES_DATA; 


long double compute_loglik_neg_binomial(double *p_model) {
    long double loglik = 0.0L;
    const long double MIN_PROB = 1e-300L; 

    for (int week = 0; week < NWEEKS_FIT; week++) {
        double prob = gsl_ran_negative_binomial_pdf(
                          CASES_DATA[week], 
                          p_model[week], 
                          OVERDISP);

	if (prob < MIN_PROB) prob = MIN_PROB;  

        loglik += logl(prob);   
	if(VERBOSE==1){
	  fprintf(stderr,"week=%d\tcases_data=%u\tp_model=%f\tloglik=%Le\n",
		week, CASES_DATA[week], p_model[week], logl(prob));
	}
    }

    return loglik;
}





void mcmc(){
  int isim;
  char * output_file;
  output_file = (char*)calloc(1000,sizeof(char)); 
  char * directory_output;
  directory_output = (char*)calloc(1000,sizeof(char)); 
  FILE *fout;

  FILE *fout2;
  double diff;
  long double a;
  //int i;

  double * p_model_inc;

  
  p_model_inc=(double *)calloc(NWEEKS_FIT, sizeof(double));
  
  sprintf(directory_output, "%soutput_files",EXP_STRING);
  fprintf(stderr,"%s\n",directory_output);
  
  if(mkdir(directory_output,S_IRWXU) != 0 && errno != EEXIST){
    perror("mkdir failed");
    exit(1);
  }
 


  sprintf(output_file, "%soutput_files/sigmas",EXP_STRING);
  fout2=fopen(output_file,"w");
  if(fout2==NULL){
    fprintf(stderr,"Error opening file sigmas for writing\n");
    exit(1);
  }

  fprintf(fout2,"%e\t%e\t%e\n",SD_BETA,SD_OVERDISP,SD_SEEDINF);
  fclose(fout2);
  

  sprintf(output_file, "%soutput_files/mcmcprogress",EXP_STRING);
  fout=fopen(output_file,"w");
  if(fout==NULL){
    fprintf(stderr,"Error opening file mcmcprogress for writing\n");
    exit(1);
  }

  

  for(isim=0;isim<NSIM;isim++){
    percent(isim,NSIM,1,stderr);
    
    if(isim==0){
      BETACURRENT=BETASTART;
      REPORTINGCURRENT=REPORTINGSTART;
      OVERDISPCURRENT=OVERDISPSTART;
      SEEDINFCURRENT=SEEDINFSTART;
  

      BETA=BETACURRENT;
      REPORTING=REPORTINGCURRENT;
      OVERDISP=OVERDISPCURRENT;
      SEEDINF=SEEDINFCURRENT;

      
      IDCURRENT=isim;
      strategy_0(isim,p_model_inc);

      
      LKHCURRENT=compute_loglik_neg_binomial(p_model_inc);
      fprintf(stderr,"initial likelihood=%Le\n",LKHCURRENT);
   
    }
    
    if(isim>0){
      BETACANDIDATE=-1.;
      REPORTINGCANDIDATE=-1.;
      OVERDISPCANDIDATE=-1.;
      SEEDINFCANDIDATE=-1.;


      while(BETACANDIDATE<0){
	BETACANDIDATE=BETACURRENT + (SD_BETA*gsl_ran_gaussian(R_GLOBAL,1.));
      }

      while(REPORTINGCANDIDATE<0 || REPORTINGCANDIDATE>1){
      	REPORTINGCANDIDATE=REPORTINGCURRENT + (SD_REPORTING*gsl_ran_gaussian(R_GLOBAL,1.));
      }
      
      while(OVERDISPCANDIDATE<0){
	OVERDISPCANDIDATE=OVERDISPCURRENT + (SD_OVERDISP*gsl_ran_gaussian(R_GLOBAL,1.));
      }
      while(SEEDINFCANDIDATE<0){
	SEEDINFCANDIDATE=SEEDINFCURRENT + (SD_SEEDINF*gsl_ran_gaussian(R_GLOBAL,1.));
      }


      
      BETA=BETACANDIDATE;
      REPORTING=REPORTINGCANDIDATE;
      OVERDISP=OVERDISPCANDIDATE;
      SEEDINF=SEEDINFCANDIDATE;

      
      strategy_0(isim,p_model_inc);

      
      LKHCANDIDATE=compute_loglik_neg_binomial(p_model_inc);
      
      // difference log-likelihood
      diff = LKHCANDIDATE - LKHCURRENT;

      // acceptance ratio
      a = expl(diff);
   
      if(gsl_ran_flat(R_GLOBAL,0,1)<=a){
	fprintf(stdout,"\taccepted!\n");
	BETACURRENT=BETACANDIDATE;
	REPORTINGCURRENT=REPORTINGCANDIDATE;
	OVERDISPCURRENT=OVERDISPCANDIDATE;
	SEEDINFCURRENT=SEEDINFCANDIDATE;
	LKHCURRENT=LKHCANDIDATE;
	IDCURRENT=isim;
      }
    }
    
    fprintf(fout,"%d\t%e\t%e\t%e\t%e\t%Le\t%d\n",isim,BETACURRENT,REPORTINGCURRENT,SEEDINFCURRENT,OVERDISPCURRENT,LKHCURRENT,IDCURRENT);
    fflush(fout);
    fflush(fout);

  }


  free(output_file);
  fclose(fout);




  free(p_model_inc);
  
  free(directory_output);

  return;
}
