#include "print.h"

char filename_EAM[128];
int error_map;

NameList nameList[] = {
	NameI (cellSize),
	NameR (threshold),
	NameR (T_max),
	NameI (NeighbourMax),
	NameI (Nelems),
	NameR (NLOOP),
	NameR (stepAvg),
	NameR (stepMovie),
	NameI (stepVacancy),
	NameI (stepAtom),
	NameI (nthreads),
	NameI (seed),
	//the max neighbour number
	NameI (nebrTabFac),
	//RDF
	NameI (limitRdf),
	NameR (rangeRdf),
	NameI (sizeHistRdf),
	NameI (stepRdf),
	NameI (error_map),
};

int GetNameList (int argc, char **argv)
{
	int id, j, k, match, ok, n;
	char buff[80], *token, filename[128], line[1024];
	FILE *fp, *input;

	sprintf (buff, "../in/run.in");
	if ((fp = fopen (buff, "r")) == 0) return (0);
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++)
		nameList[k].vStatus = 0;
	ok = 1;
	while (1){
		fgets (buff, 80, fp);
		if (feof (fp)) break;
		token = strtok (buff, "\t\n");
		if (! token) break;
		match = 0;
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k++){
			if (strcmp (token, nameList[k].vName) == 0){
				match = 1;
				if (nameList[k].vStatus == 0){
					nameList[k].vStatus = 1;
					for (j = 0; j < nameList[k].vLen; j++){
						token = strtok (NULL, ", \t\n");
						if (token){
							switch (nameList[k].vType){
								case N_I:
									*NP_I = atol (token);
									break;
								case N_R:
									*NP_R = atof (token);
									break;
							}
						} else{
							nameList[k].vStatus = 2;
							ok = 0;
						}
					}
					token = strtok (NULL, ", \t\n");
					if (token){
						nameList[k].vStatus = 3;
						ok = 0;
					}
					break;
				}else{
					nameList[k].vStatus = 4;
					ok = 0;
				}
			}
		}
		if (! match) ok = 0;
	}
	fclose (fp);
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++){
		if (nameList[k].vStatus != 1) ok = 0;
	}

	AllocMem(ith_nebrR, Nelems, double);
	AllocMem (C_threshold, Nelems, double);
	AllocMem2 (elem_preci, Nelems, 10, char);
	for (n = 0; n < Nelems; n ++) {
		ith_nebrR[n] = 0;
		C_threshold[n] = 0.;
	}

	sprintf (filename, "../in/run.in");
	input = ReadFile (filename);

	printf ("\nNameList -- string\n");
	while (1) {
		if (fgets(line, 1024, input) == NULL) break;

		GetCharVariable (line, "eam_file", filename_EAM);
		GetCharVariable (line, "MC_method", MC_method);
		GetCharVariable (line, "cell_order", cell_order);
		GetCharVariables (line, "elem_preci", elem_preci);
		GetDoubleVariables (line, "C_threshold", C_threshold);
		GetDoubleVariables (line, "ith_nebrR", ith_nebrR);

		if (strncasecmp (line, "elemName", 8) == 0) {
			strcpy (line, line+8);
			n = 0;
			while (line[0] == ' ' || line[0] == 9) strcpy (line, line+1);
			while (1) {
				sscanf (line, "%s", elemName[n]);
				while (line[0] != ' ' && line[0] != '\0') strcpy (line, line+1);
				n ++;
				while (line[0] == ' ') strcpy (line, line+1);
				if (line[0] < 65 || line[0] > 122) break;
			}
			Nelems = n;
		}
	}
	
	fclose (input);
	return (ok);
}

void PrintNameList (FILE *fp)
{
	int j, k;

	fprintf (fp, "\nNameList -- data\n");
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++){
		fprintf (fp, "%s\t", nameList[k].vName);
		if (strlen (nameList[k].vName) < 8) fprintf (fp, "\t");
		if (nameList[k].vStatus > 0){
			for (j = 0; j < nameList[k].vLen; j ++){
				switch (nameList[k].vType){
					case N_I:
						fprintf (fp, "%d ", *NP_I);
						break;
					case N_R:
						fprintf (fp, "%#g ", *NP_R);
						break;
				}
			}
		}
		switch (nameList[k].vStatus){
			case 0:
				fprintf (fp, "** no data");
				break;
			case 1:
				break;
			case 2:
				fprintf (fp, "** missing data");
				break;
			case 3:
				fprintf (fp, "** extra data");
				break;
			case 4:
				fprintf (fp, "** multiply defined");
				break;
		}
		fprintf (fp, "\n");
	}
	fprintf (fp, "----\n");
}

FILE *ReadFile (char filename[128])
{
	FILE *input;

	if ((input = fopen (filename, "r")) == NULL) {
		printf ("%s doesn't exsit\n", filename);
		exit (1);
	}

	return (input);
}

FILE *WriteFile (char filename[128])
{
	FILE *output;

	if ((output = fopen (filename, "a+")) == NULL) {
		printf ("creating %s failed\n", filename);
		exit (1);
	}

	return (output);
}

void GetCharVariable (char *line, char *name, char *variable)
{
	int size;

	if (variable[0] != 0) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);
		sscanf (line, " %s", variable);
	}

	if (variable[0] != 0) printf ("%s %s\n", name, variable);
}

int GetCharVariables (char *line, char *name, char **variable)
{
	int size, i, N;
	char *words;

	if (variable[0][0] >= 65 && variable[0][0] <= 122) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL) break;
			strcpy (variable[N], words);
			N ++;
		}
	}
	if (variable[0][0] >= 65 && variable[0][0] <= 122) {
		printf ("%s", name);
		for (i = 0; i < N; i ++) printf (" %s", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}

int GetDoubleVariables (char *line, char *name, double *variable) //the initial value of variable must be zero
{
	int size, i, N;
	char *words;

	if (fabs(variable[0]) > 1.e-60) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL) break;
			variable[N] = atof (words);
			N ++;
		}
	}
	if (fabs(variable[0]) > 1.e-60) {
		printf ("%s", name);
		for (i = 0; i < N; i ++) printf (" %.3f", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}

int GetIntVariables (char *line, char *name, int *variable) //the initial value of variable must be zero
{
	int size, i, N;
	char *words;

	if (fabs(variable[0]) != 0) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL) break;
			variable[N] = atoi (words);
			N ++;
		}
	}
	if (fabs(variable[0]) != 0) {
		printf ("%s", name);
		for (i = 0; i < N; i ++) printf (" %d", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}
