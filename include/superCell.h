#ifndef SUPERCELL_H
#define SUPERCELL_H

#include "main.h"

typedef struct {
	int Natoms;
	VecR *r;
	char name[8];
} Elem;

extern VecI cellSize;
extern char elemName[20][8];

void CreatSuperCell ();

#endif
