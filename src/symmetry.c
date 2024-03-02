#include "symmetry.h"

char pointGroup[128];
char symbol[48][128];
double ***mat;

void ObtainSym ()
{
	int i, j, k;

	if (strcmp (pointGroup,"1") == 0) {
		Nsym = 1;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
	}else if (strcmp (pointGroup, "222") == 0) {
		Nsym = 4;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 1);
		MatCopy(sym, mat, 2, 2);
		MatCopy(sym, mat, 3, 3);
	}else if (strcmp (pointGroup, "mm2") == 0) {
		Nsym = 4;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 1);
		MatCopy(sym, mat, 2, 26);
		MatCopy(sym, mat, 3, 27);
	}else if (strcmp (pointGroup, "mmm") == 0) {
		Nsym = 8;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 1);
		MatCopy(sym, mat, 2, 2);
		MatCopy(sym, mat, 3, 3);
		MatCopy(sym, mat, 4, 24);
		MatCopy(sym, mat, 5, 25);
		MatCopy(sym, mat, 6, 26);
		MatCopy(sym, mat, 7, 27);
	}else if (strcmp (pointGroup, "4") == 0) {
		Nsym = 4;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 18);
		MatCopy(sym, mat, 2, 21);
		MatCopy(sym, mat, 3, 1);
	}else if (strcmp (pointGroup, "-4") == 0) {
		Nsym = 4;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 24);
		MatCopy(sym, mat, 1, 42);
		MatCopy(sym, mat, 2, 45);
		MatCopy(sym, mat, 3, 1);
	}else if (strcmp (pointGroup, "4/m") == 0) {
		Nsym = 8;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 18);
		MatCopy(sym, mat, 2, 21);
		MatCopy(sym, mat, 3, 1);
		MatCopy(sym, mat, 4, 24);
		MatCopy(sym, mat, 5, 42);
		MatCopy(sym, mat, 6, 45);
		MatCopy(sym, mat, 7, 25);
	}else if (strcmp (pointGroup, "422") == 0) {
		Nsym = 8;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 18);
		MatCopy(sym, mat, 2, 21);
		MatCopy(sym, mat, 3, 1);
		MatCopy(sym, mat, 4, 2);
		MatCopy(sym, mat, 5, 3);
		MatCopy(sym, mat, 6, 12);
		MatCopy(sym, mat, 7, 15);
	}else if (strcmp (pointGroup, "4mm") == 0) {
		Nsym = 8;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 18);
		MatCopy(sym, mat, 2, 21);
		MatCopy(sym, mat, 3, 1);
		MatCopy(sym, mat, 4, 25);
		MatCopy(sym, mat, 5, 27);
		MatCopy(sym, mat, 6, 36);
		MatCopy(sym, mat, 7, 38);
	}else if (strcmp (pointGroup, "-42m") == 0) {
		Nsym = 8;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 24);
		MatCopy(sym, mat, 1, 42);
		MatCopy(sym, mat, 2, 45);
		MatCopy(sym, mat, 3, 1);
		MatCopy(sym, mat, 4, 2);
		MatCopy(sym, mat, 5, 3);
		MatCopy(sym, mat, 6, 36);
		MatCopy(sym, mat, 7, 38);
	}else if (strcmp (pointGroup, "4/mmm") == 0) {
		Nsym = 16;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 18);
		MatCopy(sym, mat, 2, 21);
		MatCopy(sym, mat, 3, 1);
		MatCopy(sym, mat, 4, 2);
		MatCopy(sym, mat, 5, 3);
		MatCopy(sym, mat, 6, 12);
		MatCopy(sym, mat, 7, 14);
		MatCopy(sym, mat, 8, 24);
		MatCopy(sym, mat, 9, 42);
		MatCopy(sym, mat, 10, 45);
		MatCopy(sym, mat, 11, 25);
		MatCopy(sym, mat, 12, 26);
		MatCopy(sym, mat, 13, 27);
		MatCopy(sym, mat, 14, 36);
		MatCopy(sym, mat, 15, 38);
	}else if (strcmp (pointGroup, "23") == 0) {
		Nsym = 12;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 4);
		MatCopy(sym, mat, 2, 5);
		MatCopy(sym, mat, 3, 6);
		MatCopy(sym, mat, 4, 7);
		MatCopy(sym, mat, 5, 8);
		MatCopy(sym, mat, 6, 9);
		MatCopy(sym, mat, 7, 10);
		MatCopy(sym, mat, 8, 11);
		MatCopy(sym, mat, 9, 1);
		MatCopy(sym, mat, 10, 2);
		MatCopy(sym, mat, 11, 3);
	}else if (strcmp (pointGroup, "m3") == 0) {
		Nsym = 24;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 4);
		MatCopy(sym, mat, 2, 5);
		MatCopy(sym, mat, 3, 6);
		MatCopy(sym, mat, 4, 7);
		MatCopy(sym, mat, 5, 8);
		MatCopy(sym, mat, 6, 9);
		MatCopy(sym, mat, 7, 10);
		MatCopy(sym, mat, 8, 11);
		MatCopy(sym, mat, 9, 1);
		MatCopy(sym, mat, 10, 2);
		MatCopy(sym, mat, 11, 3);
		MatCopy(sym, mat, 12, 24);
		MatCopy(sym, mat, 13, 28);
		MatCopy(sym, mat, 14, 29);
		MatCopy(sym, mat, 15, 30);
		MatCopy(sym, mat, 16, 31);
		MatCopy(sym, mat, 17, 32);
		MatCopy(sym, mat, 18, 33);
		MatCopy(sym, mat, 19, 34);
		MatCopy(sym, mat, 20, 35);
		MatCopy(sym, mat, 21, 25);
		MatCopy(sym, mat, 22, 26);
		MatCopy(sym, mat, 23, 27);
	}else if (strcmp (pointGroup, "432") == 0) {
		Nsym = 24;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 4);
		MatCopy(sym, mat, 2, 5);
		MatCopy(sym, mat, 3, 6);
		MatCopy(sym, mat, 4, 7);
		MatCopy(sym, mat, 5, 8);
		MatCopy(sym, mat, 6, 9);
		MatCopy(sym, mat, 7, 10);
		MatCopy(sym, mat, 8, 11);
		MatCopy(sym, mat, 9, 1);
		MatCopy(sym, mat, 10, 2);
		MatCopy(sym, mat, 11, 3);
		MatCopy(sym, mat, 12, 12);
		MatCopy(sym, mat, 13, 13);
		MatCopy(sym, mat, 14, 14);
		MatCopy(sym, mat, 15, 15);
		MatCopy(sym, mat, 16, 16);
		MatCopy(sym, mat, 17, 17);
		MatCopy(sym, mat, 18, 18);
		MatCopy(sym, mat, 19, 19);
		MatCopy(sym, mat, 20, 20);
		MatCopy(sym, mat, 21, 21);
		MatCopy(sym, mat, 22, 22);
		MatCopy(sym, mat, 23, 23);
	}else if (strcmp (pointGroup, "43m") == 0) {
		Nsym = 24;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 4);
		MatCopy(sym, mat, 2, 5);
		MatCopy(sym, mat, 3, 6);
		MatCopy(sym, mat, 4, 7);
		MatCopy(sym, mat, 5, 8);
		MatCopy(sym, mat, 6, 9);
		MatCopy(sym, mat, 7, 10);
		MatCopy(sym, mat, 8, 11);
		MatCopy(sym, mat, 9, 1);
		MatCopy(sym, mat, 10, 2);
		MatCopy(sym, mat, 11, 3);
		MatCopy(sym, mat, 12, 36);
		MatCopy(sym, mat, 13, 37);
		MatCopy(sym, mat, 14, 38);
		MatCopy(sym, mat, 15, 39);
		MatCopy(sym, mat, 16, 40);
		MatCopy(sym, mat, 17, 41);
		MatCopy(sym, mat, 18, 18);
		MatCopy(sym, mat, 19, 19);
		MatCopy(sym, mat, 20, 20);
		MatCopy(sym, mat, 21, 21);
		MatCopy(sym, mat, 22, 22);
		MatCopy(sym, mat, 23, 23);
	}else if (strcmp (pointGroup, "m3m") == 0) {
		Nsym = 48;
		AllocMem3 (sym, Nsym, 3, 3, double);
		MatCopy(sym, mat, 0, 0);
		MatCopy(sym, mat, 1, 4);
		MatCopy(sym, mat, 2, 5);
		MatCopy(sym, mat, 3, 6);
		MatCopy(sym, mat, 4, 7);
		MatCopy(sym, mat, 5, 8);
		MatCopy(sym, mat, 6, 9);
		MatCopy(sym, mat, 7, 10);
		MatCopy(sym, mat, 8, 11);
		MatCopy(sym, mat, 9, 1);
		MatCopy(sym, mat, 10, 2);
		MatCopy(sym, mat, 11, 3);
		MatCopy(sym, mat, 12, 12);
		MatCopy(sym, mat, 13, 13);
		MatCopy(sym, mat, 14, 14);
		MatCopy(sym, mat, 15, 15);
		MatCopy(sym, mat, 16, 16);
		MatCopy(sym, mat, 17, 17);
		MatCopy(sym, mat, 18, 18);
		MatCopy(sym, mat, 19, 19);
		MatCopy(sym, mat, 20, 20);
		MatCopy(sym, mat, 21, 21);
		MatCopy(sym, mat, 22, 22);
		MatCopy(sym, mat, 23, 23);
		MatCopy(sym, mat, 24, 24);
		MatCopy(sym, mat, 25, 28);
		MatCopy(sym, mat, 26, 29);
		MatCopy(sym, mat, 27, 30);
		MatCopy(sym, mat, 28, 31);
		MatCopy(sym, mat, 29, 32);
		MatCopy(sym, mat, 30, 33);
		MatCopy(sym, mat, 31, 34);
		MatCopy(sym, mat, 32, 35);
		MatCopy(sym, mat, 33, 25);
		MatCopy(sym, mat, 34, 26);
		MatCopy(sym, mat, 35, 27);
		MatCopy(sym, mat, 36, 36);
		MatCopy(sym, mat, 37, 37);
		MatCopy(sym, mat, 38, 38);
		MatCopy(sym, mat, 39, 39);
		MatCopy(sym, mat, 40, 40);
		MatCopy(sym, mat, 41, 41);
		MatCopy(sym, mat, 42, 42);
		MatCopy(sym, mat, 43, 43);
		MatCopy(sym, mat, 44, 44);
		MatCopy(sym, mat, 45, 45);
		MatCopy(sym, mat, 46, 46);
		MatCopy(sym, mat, 47, 47);
	}
}

void ReadMat ()
{
	FILE *input;
	char filename[128], line[1024];
	int i, j, k;

	AllocMem3 (mat, 48, 3, 3, double);

	sprintf (filename, "../in/symmetry.in");
	input = ReadFile (filename);

	for (i = 0; i < 48; i ++) {
		fgets(symbol[i], 128, input);
//		puts (symbol[i]);

		for (j = 0; j < 3; j ++) {
			fgets(line, 1024, input);
			sscanf (line, "%lg %lg %lg", &mat[i][j][0], &mat[i][j][1], &mat[i][j][2]);
//			printf ("%.0f %.0f %.0f\n", mat[i][j][0], mat[i][j][1], mat[i][j][2]);
		}
		fgets(line, 1024, input);
//		puts(line);
	}

	fclose (input);

}

void MatCopy (double ***A, double ***B, int m, int n)
{
	int i, j;

	for (i = 0; i < 3; i ++) {
		for (j = 0; j < 3; j ++) {
			A[m][i][j] = B[n][i][j];
		}
	}
	printf ("(%d): %s", m, symbol[n]);
	for (i = 0; i < 3; i ++) {
		for (j = 0; j < 3; j ++) {
			printf ("%.0f\t", A[m][i][j]);
		}
		printf ("\n");
	}
	printf ("------------------\n");
}
