#ifndef NEBRLIST_H
#define NEBRLIST_H

#include "main.h"

#define N_OFFSET 14  
#define OFFSET_VALS						\
	{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1},	\
	{1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, 		\
	{-1,-1,1}, {0,-1,1}, {1,-1,1}}
#define VCellWrap(t)						\
	if (m2v.t >= cells.t){					\
		m2v.t = 0;					\
		shift.t = region.t;				\
	}else if (m2v.t < 0){					\
		m2v.t = cells.t - 1;				\
		shift.t = - region.t;				\
	}
#define VCellWrapAll()						\
	{VCellWrap (x);						\
	 VCellWrap (y);						\
	 VCellWrap (z);}
#define DO_CELL(j, m)						\
	for (j = cellList[m]; j >= 0; j = cellList[j])
#define VWrap(v, t)					\
	if(v.t >= 0.5 * region.t) v.t -= region.t;	\
	else if (v.t < -0.5 * region.t) v.t += region.t
#define VWrapAll(v)					\
	{VWrap (v, x);					\
	 VWrap (v, y);					\
	 VWrap (v, z);}

extern VecI cells;
extern int *cellList;
extern int nebrListFlag;
extern double dispHi, rNebrShell, rCut;
extern int nMol, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;
extern Mol *mol;

void BuildNebrList ();
void CheckNebrList ();
void GetNebrs (double rNebr);

#endif
