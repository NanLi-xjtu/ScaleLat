#ifndef MAPPING_H
#define MAPPING_H

#include "main.h"

#define PeriodBoundaryGoal(v) \
	{if (v.x >= region_goal.x / 2.) v.x -= region_goal.x; \
	else if (v.x < -region_goal.x / 2.) v.x += region_goal.x;\
	if (v.y >= region_goal.y / 2.) v.y -= region_goal.y; \
	else if (v.y < -region_goal.y / 2.) v.y += region_goal.y;\
	if (v.z >= region_goal.z / 2.) v.z -= region_goal.z; \
	else if (v.z < -region_goal.z / 2.) v.z += region_goal.z;}

#define PeriodBoundaryCE(v) \
	{if (v.x >= regionCE.x - 0.01) v.x -= regionCE.x; \
	else if (v.x < -0.01) v.x += regionCE.x;\
	if (v.y >= regionCE.y - 0.01) v.y -= regionCE.y; \
	else if (v.y < -0.01) v.y += regionCE.y;\
	if (v.z >= regionCE.z - 0.01) v.z -= regionCE.z; \
	else if (v.z < -0.01) v.z += regionCE.z;}

#define PeriodBoundaryCEdr(v) \
	{if (v.x >= regionCE.x / 2.) v.x -= regionCE.x; \
	else if (v.x < -regionCE.x / 2.) v.x += regionCE.x;\
	if (v.y >= regionCE.y / 2.) v.y -= regionCE.y; \
	else if (v.y < -regionCE.y / 2.) v.y += regionCE.y;\
	if (v.z >= regionCE.z / 2.) v.z -= regionCE.z; \
	else if (v.z < -regionCE.z / 2.) v.z += regionCE.z;}

#define VCellWrapGoal(t)						\
	if (m2v.t >= cells.t){					\
		m2v.t = 0;					\
		shift.t = region_goal.t;				\
	}else if (m2v.t < 0){					\
		m2v.t = cells.t - 1;				\
		shift.t = - region_goal.t;				\
	}
#define VCellWrapAllGoal()						\
	{VCellWrapGoal (x);						\
	 VCellWrapGoal (y);						\
	 VCellWrapGoal (z);}

typedef struct{
	int small_rc_id, big_rc_id;
	double Csmall_rc, Cbig_rc;
} MATCH;

extern int NcluMax;
extern VecI clusterSize;
extern int Nsym;
extern double ***sym;

void Mapping ();
int ClusterExtract (int N, Mol *cluster, CTABLE *C, Mol *rc, int mode);
int Symmetry (VecR *Vdr1, char **element1, VecR *Vdr2, char **element2);
void GoalClusterExtract ();
void MapStructure ();
void MapConcentrationTable ();
void MapClusterExtract ();
void PrintCluster (int N, VecR *c, char **element);
void PrintRepresentiveCluster ();
void PrintConcentration ();
void PrintMapStructure ();
void PrintMatch ();

#endif
