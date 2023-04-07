#ifndef PRINT_H
#define PRINT_H

#include "main.h"

typedef enum {N_I, N_R} VType;

//how to use NameList
/****************************
NameList nameList[] = {
	NameI (intVariable),
	NameR (realVector),
};
*****************************/
typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;

#define NameI(x) {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x) {#x, &x, N_R, sizeof (x) / sizeof (real)}
#define NP_I ((int *) (nameList[k].vPtr) + j)
#define NP_R ((real *) (nameList[k].vPtr) + j)

extern char filename_EAM[128];
extern int error_map;

int GetNameList (int argc, char **argv);
void PrintNameList (FILE *fp);
FILE *ReadFile (char filename[128]);
FILE *WriteFile (char filename[128]);
void GetCharVariable (char *line, char *name, char *variable);
int GetCharVariables (char *line, char *name, char **variable);
int GetDoubleVariables (char *line, char *name, double *variable);
int GetIntVariables (char *line, char *name, int *variable);

#endif
