#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "main.h"

#define PeriodBoundaryGoal(v) \
	if (v.x >= region_goal.x / 2.) v.x -= region_goal.x; \
	else if (v.x < -region_goal.x / 2.) v.x += region_goal.x;\
	if (v.y >= region_goal.y / 2.) v.y -= region_goal.y; \
	else if (v.y < -region_goal.y / 2.) v.y += region_goal.y;\
	if (v.z >= region_goal.z / 2.) v.z -= region_goal.z; \
	else if (v.z < -region_goal.z / 2.) v.z += region_goal.z;\

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
//RDF
extern double **histRdf, rangeRdf;
extern int countRdf, limitRdf, sizeHistRdf, stepRdf;
extern double *C_threshold, *ith_nebrR;
extern char **elem_preci;

double ShortRangeOrder (char *A, char *B, char *filename);
void EvalRdf ();
void IdenPrecipitates ();
void CalSROMovie (char *file, char *element1, char *element2);
void CalSROMovie_TOTAL ();
void CalAutocorrelationFunction ();
void ClusterExtract ();
int Symmetry (Mol mol1, Mol mol2);
void KMC_Cluster ();

#endif
