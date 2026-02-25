#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "header.h"

void save_vector_of_double(const char *filename,const double *vector,unsigned maxlen){
    FILE *fp = NULL;
    register unsigned int i;
    if ((fp = fopen(filename, "w")) == NULL) {
      fprintf(stderr, "Error opening <%s> for writing\n", filename);
      return;
    }
    for (i = 0; i < (maxlen); i++) {
      fprintf(fp, "%f\n", vector[i]);
    }
    fprintf(fp, "\n");
    fflush(fp);
    fflush(fp);

    fclose(fp);
}

void save_vector_of_int(const char *filename,const int *vector,unsigned maxlen){
    FILE *fp = NULL;
    register unsigned int i;
    if ((fp = fopen(filename, "w")) == NULL) {
      fprintf(stderr, "Error opening <%s> for writing\n", filename);
      return;
    }
    for (i = 0; i < (maxlen); i++) {
      fprintf(fp, "%d\n", vector[i]);
    }
    fprintf(fp, "\n");
    fflush(fp);
    fflush(fp);

    fclose(fp);
}

void save_matrix_of_double(const char *filename, double **matrix,unsigned nrow, unsigned ncol){
  FILE *fp = NULL;
  register unsigned int i,j;
  if ((fp = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "Error opening <%s> for writing\n", filename);
    return;
  }
  for (i = 0; i < (nrow); i++) {
    for (j = 0; j < (ncol); j++) {
      fprintf(fp, "%f\t", matrix[i][j]);
    }
    fprintf(fp, "\n");
  }

  fflush(fp);
  fflush(fp);

  fclose(fp);
}
