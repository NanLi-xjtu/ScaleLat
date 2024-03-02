#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "main.h"

extern char pointGroup[128];

void ObtainSym ();
void ReadMat ();
void MatCopy (double ***A, double ***B, int m, int n);

#endif
