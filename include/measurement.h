#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "main.h"

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

#endif
