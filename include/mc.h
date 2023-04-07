#ifndef MC_H
#define MC_H

#include "main.h"

typedef struct {
	int Nnebrs, *nebr;
	double distance;
} NebrLists;

typedef struct {
	char element[128];
	NebrLists *nebrLists;
} NebrElem;

typedef struct {
	VecR r;
	NebrLists *nebrLists;
//	NebrElem *nebrElem;
	char element[128];
	int type;
} Atom;

#define DO_ATOM for (n = 0; n < Natoms[0]; n ++)
#define DO_MOL for (n = 0; n < nMol; n ++)
#define PeriodBoundary(v) \
	if (v.x >= region.x / 2.) v.x -= region.x; \
	else if (v.x < -region.x / 2.) v.x += region.x;\
	if (v.y >= region.y / 2.) v.y -= region.y; \
	else if (v.y < -region.y / 2.) v.y += region.y;\
	if (v.z >= region.z / 2.) v.z -= region.z; \
	else if (v.z < -region.z / 2.) v.z += region.z;\


extern Atom *atom;
extern VecI cellSize;
extern VecR region;
extern double TUnit, EUnit, VUnit, lUnit, eUnit, T_max, T_min, dT, dT_min, dT_max, threshold, *pairDists, NLOOP, stepMovie, stepAvg;
extern int *Natoms, Nelems, NeighbourMax, Ndists;
extern long long Nloop, stepCount;
extern int nthreads, seed;
extern char MC_method[128], cell_order[128];
extern int stepVacancy, stepAtom;

void Setup ();
void SetupMD ();
void Restart ();
void MonteCarloVacancy (double T);
void MonteCarloExchangeAtoms (double T);
void MonteCarloMixed (double T);
void MonteCarloDoubleT (double temperature);
void PrintAtoms ();
void Get_pair_distance ();
void Get_neighbour ();
double Get_dE (int vacancy, int nebr);
double Get_dE_new (int vacancy, int nebr);
double Get_dE_OneAtomChange (int id, char *element);
double Get_dE_AtomsExchange (int m1, int m2);
void PrintNeighbour ();
void PrintMovie ();

#endif
