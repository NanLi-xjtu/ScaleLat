#ifndef DEFINE_H
#define DEFINE_H

#include "main.h"

typedef double real;
typedef struct {
	real x,y,z;
} VecR;
typedef struct {
	int x,y,z;
} VecI;
typedef struct{
	VecR r, rv, ra;
	char elem[8]; //the type of element of the atom
	int Nnebr, *nebr_id, nebr[2000];
	int preci_flag;
	int interface_flag;
	int NCEatom, *CEatom; //atoms of each cluster
	double C;
	int *Nelem;
	double *Celem;
	int id;
	int Nsym; //the number of symmetry cluters
	int match_flag;
} Mol;
typedef struct{
	double C;
	int *Nc;
} CTABLE;

#define NDIM 3
#define Nthreads 22

//constant
#define eleChar 1.60218e-19 //--C
#define eleMass 9.10956e-31 //--kg
#define kB 1.3806505e-23 //--J/K
#define epsilon0 8.854187817e-12 //--F/m
#define planck 6.62607e-34 //--Js
#define planck_ 1.05457e-34 //--Js

//the allocations of one- and two- and three- dimensional arrays of any kind of variable or structure
#define AllocMem(a, n, t) 				\
	a = (t *) malloc ((n) * sizeof (t))
#define AllocMem2(a, n1, n2, t)				\
	AllocMem (a, n1, t *);				\
	AllocMem (a[0], (n1) * (n2), t);		\
	for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;
#define AllocMem3(a, n1, n2, n3, t)					\
	AllocMem2(a, n1, n2, t *);					\
	AllocMem(a[0][0], (n1) * (n2) * (n3), t);			\
	for (k = 1; k < n1; k ++) a[k][0] = a[k - 1][0] + n2 * n3;	\
	for (j = 0; j < n1; j++){					\
		for (k = 1; k < n2; k++) a[j][k] = a[j][k - 1] + n3;	\
	}

#define ReAllocMem(a, n, t) 				\
	a = (t *) realloc (a, (n) * sizeof (t))

#define Sqr(x)	((x)*(x))
#define Cube(x)	((x)*(x)*(x))
#define Min(x1, x2)			\
	(((x1) < (x2)) ? (x1) : (x2))
#define Max(x1, x2)			\
	(((x1) > (x2)) ? (x1) : (x2))
#define MIN(x1, x2)			\
	(((x1) < (x2)) ? (x1) : (x2))
#define MAX(x1, x2)			\
	(((x1) > (x2)) ? (x1) : (x2))
#define VCMax(v)			\
	(Max(Max((v).x, (v).y), (v).z))
#define VCSum(v)			\
	((v).x + (v).y + (v).z)
#define VProd(v)	((v).x * (v).y * (v).z)
#define VDot(v1, v2)			\
	((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)
#define VLen(v)		sqrt (VDot (v, v))
#define VLenSq(v)	VDot(v, v)
#define VLinear(p, s)			\
	(((p).z * (s).y + (p).y) * (s).x + (p).x)
#define VSet(v, sx, sy, sz)		\
	(v).x = sx,			\
	(v).y = sy,			\
	(v).z = sz
#define VSetAll(v, s)	VSet(v, s, s, s)
#define VZero(v)	VSetAll(v, 0)
#define VCopy(v1,v2)			\
	(v1).x = (v2).x,		\
	(v1).y = (v2).y,		\
	(v1).z = (v2).z
#define VScale(v,s)			\
	(v).x *= s,			\
	(v).y *= s,			\
	(v).z *= s
#define VSCopy(v2, s1, v1)		\
	(v2).x = (s1) * (v1).x,		\
	(v2).y = (s1) * (v1).y,		\
	(v2).z = (s1) * (v1).z
#define VAdd(v1, v2, v3)		\
	(v1).x = (v2).x + (v3).x,	\
	(v1).y = (v2).y + (v3).y,	\
	(v1).z = (v2).z + (v3).z
#define VSAdd(v1, v2, s3, v3)		\
	(v1).x = (v2).x + (s3) * (v3).x,\
	(v1).y = (v2).y + (s3) * (v3).y,\
	(v1).z = (v2).z + (s3) * (v3).z
#define VVAdd(v1, v2)		VAdd(v1, v1, v2)
#define VVSAdd(v1, s2, v2)	VSAdd(v1, v1, s2, v2)
#define VSub(v1, v2, v3)		\
	(v1).x = (v2).x - (v3).x,	\
	(v1).y = (v2).y - (v3).y,	\
	(v1).z = (v2).z - (v3).z	
#define VVSub(v1, v2)	VSub (v1, v1, v2)
#define VMul(v1, v2, v3)		\
	(v1).x = (v2).x * (v3).x,	\
	(v1).y = (v2).y * (v3).y,	\
	(v1).z = (v2).z * (v3).z	
#define VDiv(v1, v2, v3)		\
	(v1).x = (v2).x / (v3).x,	\
	(v1).y = (v2).y / (v3).y,	\
	(v1).z = (v2).z / (v3).z

#endif
