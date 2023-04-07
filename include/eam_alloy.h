#ifndef EAM_ALLOY_H
#define EAM_ALLOY_H

#include "main.h"

typedef struct {
  char **elements;
  int nelements,nrho,nr;
  double drho,dr,cut;
  double *mass;
  double **frho,**rhor,***z2r;
} Setfl;

extern Setfl *setfl;

void coeff_alloy(char *filename);
void read_file_alloy(char *filename);
void file2array_alloy ();
void init_style_alloy();
void printfeam_alloy ();
void printfspline_alloy (double ***xxxspline);
void ComputeForcesEamPoten_alloy ();

#endif
