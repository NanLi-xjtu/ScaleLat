#include "superCell.h"

VecI cellSize;
char elemName[20][8];

float Mould(double x, double y, double z)
{
	return sqrt(x*x+y*y+z*z);
}   

void CreatSuperCell ()
{
	char elemName_CONTCAR[20][8];
	Elem *unce, *suce;
	FILE *input, *output;
	char line[1024], cellName[128], filename[128], command[128];
	int i, j, k, l, m, n, row, sum, flag;
	int NatomsUnCe, NatomsSuCe;
	double scaling, xdata, ydata, zdata;
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, a, b, c; //lattice parameters
	float X_1, Y_1, Z_1; // Atomic positions in Cartesian coordinates
	float G_1x, G_1y, G_1z; // Atomic positions in supercell
	VecR region;

	sprintf (command, "mkdir ../out/superCell/%d_%d_%d", cellSize.x, cellSize.y, cellSize.z);
	system (command);

	//read input file
	input = ReadFile ("../in/CONTCAR");
	row = 1;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;
		switch (row) {
			case 1: 
				sscanf (line, "%s", cellName); 
				puts (cellName);
				break;
			case 2: 
				sscanf (line, "%lg", &scaling); 
				printf ("%f\n", scaling);
				break;
			case 3: 
				sscanf (line, "%lg %lg %lg", &a1, &a2, &a3);
				printf ("%f\t%f\t%f\n", a1, a2, a3);
				break;
			case 4: 
				sscanf (line, "%lg %lg %lg", &b1, &b2, &b3);
				printf ("%f\t%f\t%f\n", b1, b2, b3);
				break;
			case 5: 
				sscanf (line, "%lg %lg %lg", &c1, &c2, &c3);
				printf ("%f\t%f\t%f\n", c1, c2, c3);
				break;
			case 6:
				n = 0;
				while (line[0] == ' ') strcpy (line, line+1);
				while (1) {
					sscanf (line, "%s", elemName_CONTCAR[n]);
					while (line[0] != ' ' && line[0] != '\0') strcpy (line, line+1);
					n ++;
					while (line[0] == ' ') strcpy (line, line+1);
					if (line[0] < 65 || line[0] > 122) break;
				}
/*
				if (Nelems != n) {
					printf ("error(superCell.c): Nelems\n");
					exit (1);
				}
*/
				for (n = 0; n < Nelems; n ++) {
					printf("%s ", elemName_CONTCAR[n]);
					flag = 0;
					for (i = 0; i < Nelems; i ++) {
						if (strcmp(elemName_CONTCAR[n], elemName[i]) == 0) flag = 1;
					}
/*
					if (flag == 0) {
						printf ("error(superCell.c): elemName\n");
						exit (1);
					}
*/
				}
				printf ("\n");
				AllocMem (unce, Nelems, Elem);
				AllocMem (suce, Nelems, Elem);
				break;
			case 7:
				for (n = 0; n < Nelems; n ++) {
					while (line[0] == ' ') strcpy (line, line+1);
					sscanf (line, "%d", &unce[n].Natoms);
					while (line[0] != ' ' && line[0] != '\0') strcpy (line, line+1);
					if (line[0] == '\0') break;
				}
				for (n = 0; n < Nelems; n ++) {
					AllocMem (unce[n].r, unce[n].Natoms, VecR);
					printf("%d ", unce[n].Natoms);
				}
				printf ("\n");
				break;
			case 8: i = j = k = 0; break;
			default: 
				if (i == unce[j].Natoms) {
					j ++;
					k = i = 0;
				}
				sscanf (line, "%lg %lg %lg", &unce[j].r[k].x, &unce[j].r[k].y, &unce[j].r[k].z);
				i ++;
				k ++;
				break;
		}
		row ++;
	}
	fclose (input);

	//output unit cell
	sprintf(filename,"../out/superCell/%d_%d_%d/Unit_cell_coordinates", cellSize.x, cellSize.y, cellSize.z);
	output = WriteFile (filename);
	for (i = 0; i < Nelems; i ++) {
		for (j = 0; j < unce[i].Natoms; j ++) {
			xdata = unce[i].r[j].x;
			ydata = unce[i].r[j].y;
			zdata = unce[i].r[j].z;
			X_1 = xdata * a1 + ydata * b1 + zdata * c1;
			Y_1 = xdata * a2 + ydata * b2 + zdata * c2;
			Z_1 = xdata * a3 + ydata * b3 + zdata * c3;
			printf ("%f\t%f\t%f\n", xdata, ydata, zdata);
			fprintf (output, "%f\t%f\t%f\n", X_1, Y_1, Z_1);
		}
	}
	fclose (output);

	//super cell initialize
	printf("The size of supercell are defined by three integers:\n");
	printf("Enter the size of supercell:%d %d %d\n", cellSize.x, cellSize.y, cellSize.z);
	sum = 0;
	for (i = 0; i < Nelems; i ++) sum += unce[i].Natoms;
	NatomsUnCe = sum; //the number of unit cell atoms
	NatomsSuCe = NatomsUnCe * VProd (cellSize);
	for (i = 0; i < Nelems; i ++) {
		suce[i].Natoms = unce[i].Natoms * VProd (cellSize);
		AllocMem (suce[i].r, suce[i].Natoms, VecR);
	}
	//creat super cell
	for (m = 0; m < Nelems; m ++) {
		i = 0;
		for (n = 0; n < unce[m].Natoms; n ++) {
			for(j=0;j<cellSize.x;j++){
				for(k=0;k<cellSize.y;k++){
					for(l=0;l<cellSize.z;l++){
						xdata = unce[m].r[n].x + j;
						ydata = unce[m].r[n].y + k;
						zdata = unce[m].r[n].z + l;
						X_1 = xdata * a1 + ydata * b1 + zdata * c1;
						Y_1 = xdata * a2 + ydata * b2 + zdata * c2;
						Z_1 = xdata * a3 + ydata * b3 + zdata * c3;
						VSet (suce[m].r[i], X_1, Y_1, Z_1);
						VScale (suce[m].r[i], scaling);
						i ++;
					}
				}
			}
		}
	}

	//output supercell lattice
	a = Mould (a1, a2, a3);
	b = Mould (b1, b2, b3);
	c = Mould (c1, c2, c3);
	sprintf(filename,"../out/superCell/%d_%d_%d/Lattice_supercell", cellSize.x, cellSize.y, cellSize.z);
	output = WriteFile (filename);
	region.x = cellSize.x * a;
	region.y = cellSize.y * b;
	region.z = cellSize.z * c;
	fprintf (output, "%f\n%f\n%f\n", region.x, region.y, region.z);
	fclose (output);

	//output super cell
	sprintf(filename,"../out/superCell/%d_%d_%d/Supercell.xyz", cellSize.x, cellSize.y, cellSize.z);
	output = WriteFile (filename);
	fprintf(output,"%d\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" SolutionReader properties=species:S:1:pos:R:3\n", NatomsSuCe, region.x, region.y, region.z);
	k = 0;
	for (i = 0; i < Nelems; i ++) {
		for (j = 0; j < suce[i].Natoms; j ++) {
			fprintf (output, "%s %f %f %f\n", elemName_CONTCAR[i], suce[i].r[j].x, suce[i].r[j].y, suce[i].r[j].z);
			k ++;
		}
	}
	fclose (output);
	sprintf (command, "cp ../out/superCell/%d_%d_%d/Supercell.xyz ../out/superCell/Supercell.xyz", cellSize.x, cellSize.y, cellSize.z);
	system (command);

	//output POSCAR
	sprintf(filename,"../out/superCell/%d_%d_%d/POSCAR", cellSize.x, cellSize.y, cellSize.z);
	output = WriteFile (filename);
	fprintf (output, "%s\n", cellName);
	fprintf (output, "\t1.0\n"); 
	fprintf (output, "\t%f\t%f\t%f\n", scaling*a1*cellSize.x, scaling*a2*cellSize.x, scaling*a3*cellSize.x);
	fprintf (output, "\t%f\t%f\t%f\n", scaling*b1*cellSize.y, scaling*b2*cellSize.y, scaling*b3*cellSize.y);
	fprintf (output, "\t%f\t%f\t%f\n", scaling*c1*cellSize.z, scaling*c2*cellSize.z, scaling*c3*cellSize.z);
	for (i = 0; i < Nelems; i ++) fprintf (output, "%s ", elemName[i]);
	fprintf (output, "\n");
	for (i = 0; i < Nelems; i ++) fprintf (output, "%d ", suce[i].Natoms);
	fprintf(output,"\nDirect\n");
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < unce[m].Natoms; n ++) {
			for(j=0;j<cellSize.x;j++){
				for(k=0;k<cellSize.y;k++){
					for(l=0;l<cellSize.z;l++){
						xdata = unce[m].r[n].x + j;
						ydata = unce[m].r[n].y + k;
						zdata = unce[m].r[n].z + l;
						fprintf (output, "     %f     %f     %f\n", \
							 xdata/cellSize.x, ydata/cellSize.y, zdata/cellSize.z);
					}
				}
			}
		}
	}
	fclose (output);

	free (unce);
	free (suce);
	
}
