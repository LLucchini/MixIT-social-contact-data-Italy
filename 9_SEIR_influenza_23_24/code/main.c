#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include "header.h"
#include "global.h"
#include <sys/stat.h>


int main(int argc, char *argv[]){ 
  /*read command line*/
  char * char_seed;
  char * char_fit;
  char **flag,**value;
  int nflags;
  int exit_parser; 
  int i;
  char *param_file;
  char *char_directory=NULL;
  char help[]="\nhelp FRvar:\n  -s [seed:0] -rerun [rerun:0] -verbose [verbose:1] -e [name start directory]";

  char_directory=(char*)calloc(1000,sizeof(char));
  param_file=(char*)calloc(1000,sizeof(char));
  
  /*read command line parameters*/
  if(argc<2){
    fprintf(stderr,"%s\n",help);
    exit(1);
  }
  if(argc==2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-help")==0)){
    fprintf(stderr,"%s\n",help);
    exit(1);
  }
  if((exit_parser=parser(argc,argv,&flag,&value,&nflags))!=0){
    fprintf(stderr,"parser: syntax error\n");
    fprintf(stderr,"%s\n",help);
    exit(1);
  }
  if((char_seed=get_value(flag,value,nflags,"-s"))==NULL)
    char_seed="0";

  
  ISEED=atoi(char_seed);

  EXP_STRING = (char*)calloc(1000,sizeof(char)); 


  if((char_directory=get_value(flag,value,nflags,"-e"))== NULL)
    char_directory="";
  else
    strcat(char_directory, "/");
  

  if((char_fit=get_value(flag,value,nflags,"-verbose"))== NULL)
    VERBOSE=1;
  else
    VERBOSE=atoi(char_fit);
  fprintf(stderr, "VERBOSE=%d \n", VERBOSE);

  
  if((char_fit=get_value(flag,value,nflags,"-rerun"))== NULL)
    RERUN=0;
  else
    RERUN=atoi(char_fit);

  fprintf(stderr, "RERUN=%d \n", RERUN);


  
  if(VERBOSE==1){
    fprintf(stderr, "%s \n", char_directory);
  
    fprintf(stdout,"%s\n\n\n",char_directory);

    fprintf(stdout,"%s\n\n\n",EXP_STRING);

    fprintf(stdout,"seed=%d\n",ISEED);
  }
  
  /*initialize gsl*/
  setenv("GSL_RNG_SEED",  char_seed, 1);
  gsl_rng_env_setup();
  R_GLOBAL = gsl_rng_alloc (gsl_rng_default);

  /* dynamic allocation of memory */

  C_MATRIX_PRESENCE=(double**)calloc(NAGES,sizeof(double*));
  for(i=0;i<NAGES;i++){
    C_MATRIX_PRESENCE[i]=(double*)calloc(NAGES,sizeof(double));
  }

  C_MATRIX_NOT_PRESENCE=(double**)calloc(NAGES,sizeof(double*));
  for(i=0;i<NAGES;i++){
    C_MATRIX_NOT_PRESENCE[i]=(double*)calloc(NAGES,sizeof(double));
  }


  /*Epidemiological data*/
  DATA_INC=(double*)calloc(NWEEKS_FIT,sizeof(double));
  CASES_DATA=(unsigned int*)calloc(NWEEKS_FIT,sizeof(unsigned int));

  POP_AGE_INPUT=(double*)calloc(NAGES,sizeof(double));
  POP_AGE_INPUT_P=(double*)calloc(NAGES,sizeof(double));
  POP_AGE_INPUT_NP=(double*)calloc(NAGES,sizeof(double));
  SUSC_AGE_INPUT=(double*)calloc(NAGES,sizeof(double));

  
  sprintf(EXP_STRING, "%s", char_directory);

  if(VERBOSE==1){
    fprintf(stderr,"EXP_STRING=%s\n", EXP_STRING);
  }
  
  
  /*read parameters*/
  sprintf(param_file,"%s/input_files/parameters",EXP_STRING);
  if(read_param(param_file,'\t')!=0){
    fprintf(stderr,"error in read_param\n");
    exit(1);
  }
  
  fprintf(stderr,"OMEGA: %f\n",OMEGA);
  fprintf(stderr,"GAMMA: %f\n",GAMMA);
  fprintf(stderr,"NSIM: %d\n",NSIM);
  fprintf(stderr,"BETASTART: %f\n",BETASTART);
  fprintf(stderr,"REPORTINGSTART: %f\n",REPORTINGSTART);
  fprintf(stderr,"OVERDISPSTART: %f\n",OVERDISPSTART);
  fprintf(stderr,"SEEDINFSTART: %f\n",SEEDINFSTART);
  fprintf(stderr,"SD_BETA: %f\n",SD_BETA);
  fprintf(stderr,"SD_REPORTING: %f\n",SD_REPORTING);
  fprintf(stderr,"SD_OVERDISP: %f\n",SD_OVERDISP);
  fprintf(stderr,"SD_SEEDINF: %f\n",SD_SEEDINF);



  if(read_demo_data('\t')!=0){
    fprintf(stderr,"error in read_ demo_data\n");
    exit(1);
  }

  if(RERUN==0){
    /*Run calibration with mcmc*/
    if(read_incidence()!=0){
      fprintf(stderr,"error in read_incidence\n");
      exit(1);
    }
    mcmc();
  }else{
    /*Rerun simulations*/
    rerun();
  }
    
  free(EXP_STRING);
  free(param_file); 
  

  
  for(i=0;i<NAGES;i++){
    free(C_MATRIX_PRESENCE[i]);
    free(C_MATRIX_NOT_PRESENCE[i]);
  }
  free(C_MATRIX_PRESENCE);
  free(C_MATRIX_NOT_PRESENCE);


  free(POP_AGE_INPUT);
  free(POP_AGE_INPUT_P);
  free(POP_AGE_INPUT_NP);
  //free(REPORTING_AGE_INPUT);
  free(SUSC_AGE_INPUT);

  free(DATA_INC);
  free(CASES_DATA);

  return 0;
}

