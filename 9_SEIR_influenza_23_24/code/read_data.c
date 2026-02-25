#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"

extern double * POP_AGE_INPUT;
extern double * POP_AGE_INPUT_P;
extern double * POP_AGE_INPUT_NP;
extern double *SUSC_AGE_INPUT;


extern double POPULATIONSIZE;
extern double ** C_MATRIX_PRESENCE;
extern double ** C_MATRIX_NOT_PRESENCE;

extern int NAGES;
extern char *EXP_STRING;

extern double *DATA_INC;
extern unsigned int *CASES_DATA;

extern int NWEEKS_FIT; 
extern int WEEK_START_FIT;
extern int WEEK_END_FIT;
extern int VERBOSE;

int read_demo_data(char sep){
  char *line;
  FILE *fp;
  int i,j;
  double tmp, tmpP,tmpNP;
  char *demo_file;

  demo_file=(char*)calloc(1000,sizeof(char));

  sprintf(demo_file,"%s/input_files/pop_by_age",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }

  tmp=0.0;
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    sscanf(line,"%lf", &POP_AGE_INPUT[i]);
    tmp+=POP_AGE_INPUT[i];
  }
  POPULATIONSIZE=tmp;
  fprintf(stderr,"POPSIZE=%lf\n",POPULATIONSIZE);

  fclose(fp);

  
  sprintf(demo_file,"%s/input_files/susc_age_input",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }
  
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    sscanf(line,"%lf", &SUSC_AGE_INPUT[i]);
    if(VERBOSE==1){
      fprintf(stderr,"SUSC_AGE=%f\n",SUSC_AGE_INPUT[i]);
    }
  }
  
  fclose(fp);

  sprintf(demo_file,"%s/input_files/pop_by_age_presence",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }
  
  tmpP=0.0;
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    sscanf(line,"%lf", &POP_AGE_INPUT_P[i]);
    tmpP+=POP_AGE_INPUT_P[i];
  }
  fclose(fp);

  sprintf(demo_file,"%s/input_files/pop_by_age_not_presence",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }
  
  tmpNP=0.0;
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    sscanf(line,"%lf", &POP_AGE_INPUT_NP[i]);
    tmpNP+=POP_AGE_INPUT_NP[i];
  }
  
  fclose(fp);

  fprintf(stderr,"POPSIZE_P=%lf\t POPSIZE_NP=%lf\n",tmpP,tmpNP);

  POPULATIONSIZE=tmpP+tmpNP;
  fprintf(stderr,"POPSIZE=%lf\n",POPULATIONSIZE);
   

  
  sprintf(demo_file,"%s/input_files/matrix_IT_presence",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }

  tmp=0.0;
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    for(j=0;j<NAGES;j++){
      sscanf(line,"%lf", &C_MATRIX_PRESENCE[i][j]);
      tmp+=C_MATRIX_PRESENCE[i][j];
      line = (char *)strchr(line, sep);
      line++;
    }
  }
  fprintf(stderr,"sum_contact_matrix_presence=%f\n",tmp); 
  fclose(fp);


  sprintf(demo_file,"%s/input_files/matrix_IT_not_presence",EXP_STRING);
  if(!(fp = fopen(demo_file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",demo_file);
    return 1;
  }

  tmp=0.0;
  for(i=0;i<NAGES;i++){
    get_line(&line,fp);
    for(j=0;j<NAGES;j++){
      sscanf(line,"%lf", &C_MATRIX_NOT_PRESENCE[i][j]);
      tmp+=C_MATRIX_NOT_PRESENCE[i][j];
      line = (char *)strchr(line, sep);
      line++;
    }
  }
  fprintf(stderr,"sum_contact_matrix_not_presence=%f\n",tmp);
  fclose(fp);


  free(demo_file);


  return 0;
}

int read_incidence() {
    FILE *fp = fopen("input_files/ILI_A_2023_2024", "r");
    if (!fp) {
        perror("Errore apertura file");
        return 1;
    }

    int settimana;
    double valore;
    int count = 0;
    int capacity = 100;

    // Array dinamici
    int *weeks = malloc(capacity * sizeof(int));
    double *values = malloc(capacity * sizeof(double));

    if (!weeks || !values) {
        perror("Errore allocazione memoria");
        return 1;
    }

    while (fscanf(fp, "%d\t%lf", &settimana, &valore) == 2) {
        if (count >= capacity) {
            capacity *= 2;
            weeks = realloc(weeks, capacity * sizeof(int));
            values = realloc(values, capacity * sizeof(double));
            if (!weeks || !values) {
                perror("Errore realloc");
                return 1;
            }
        }
        weeks[count] = settimana;
        values[count] = valore;
        count++;
    }

    fclose(fp);

    
    for (int i = 0; i < NWEEKS_FIT; i++) {
      DATA_INC[i]=values[WEEK_START_FIT+i];
      
    }

    for (int i = 0; i < NWEEKS_FIT; i++) {
      CASES_DATA[i]=(unsigned int) round(DATA_INC[i]*POPULATIONSIZE/1000.);


    }
    free(weeks);
    free(values);

    return 0;
}


int readMatrix(const char *filename,
               double *beta,
               double *reporting,
               double *seedinf,
               size_t nrows) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        return -1;
    }

    size_t i = 0;
    while (i < nrows && fscanf(fp, "%lf %lf %lf",
                               &beta[i], &reporting[i], &seedinf[i]) == 3) {
        i++;
    }

    fclose(fp);

    if (i != nrows) {
        fprintf(stderr, "Warning: only %zu rows read out of %zu\n", i, nrows);
        return -1;
    }

    return 0;
}
