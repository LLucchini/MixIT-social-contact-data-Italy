#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"

extern double OMEGA;
extern double GAMMA;
extern int NSIM;
extern double BETASTART;
extern double OVERDISPSTART;
extern double SEEDINFSTART;
extern double REPORTINGSTART;
extern double SD_BETA;
extern double SD_OVERDISP;
extern double SD_SEEDINF;  
extern double SD_REPORTING;



int read_param(char file[],char sep){

  char *line;
  int out_get_line=2;
  FILE *fp;
  char *string_who, *string_value;

  string_who=(char*)calloc(100,sizeof(char));
  string_value=(char*)calloc(100,sizeof(char));

  if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"read_param: error opening file %s for reading\n",file);
    return 1;
  }

  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_param: line of file %s does not end in newline\n",file);
	break;
      case 1:
	fprintf(stderr,"read_param: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_param: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_param: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    sscanf(line,"%s", string_who);
    line = strchr(line, sep);
    if(line) {
    line++;
    sscanf(line,"%s", string_value);
    }

    sscanf(line,"%s", string_value);
     if(strcmp(string_who,"OMEGA")==0)
      OMEGA=atof(string_value);
    if(strcmp(string_who,"GAMMA")==0)
      GAMMA=atof(string_value);
    if(strcmp(string_who,"BETASTART")==0)
      BETASTART=atof(string_value);
    if(strcmp(string_who,"REPORTINGSTART")==0)
      REPORTINGSTART=atof(string_value);
    if(strcmp(string_who,"OVERDISPSTART")==0)
      OVERDISPSTART=atof(string_value);
    if(strcmp(string_who,"SEEDINFSTART")==0)
      SEEDINFSTART=atof(string_value);
    if(strcmp(string_who,"SD_BETA")==0)
      SD_BETA=atof(string_value);
    if(strcmp(string_who,"SD_REPORTING")==0)
      SD_REPORTING=atof(string_value);
    if(strcmp(string_who,"SD_OVERDISP")==0)
      SD_OVERDISP=atof(string_value);
    if(strcmp(string_who,"SD_SEEDINF")==0)
      SD_SEEDINF=atof(string_value);
    if(strcmp(string_who,"NSIM")==0)
      NSIM=atoi(string_value);
  }

  free(string_who);
  free(string_value);


  return 2;

}
