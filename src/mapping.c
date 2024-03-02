#include "mapping.h"

MATCH *mat;
int Nmat;
Mol *goal;
int Ngoal;
VecR regionCE, region_goal;
int Nrc_goal, Nrc_map; //the number of representive clusters
Mol *rc_goal, *rc_map; //representive clusters
int NC; //the number of concentations
CTABLE *C_goal, *C_map; //concentation table
VecI clusterSize;
int Nsym;
double ***sym;

void Mapping ()
{
	int i, j, k, i1, j1, k1, i2, j2, k2, n, id, flag_elem, flag_sym, A, B, loop, NsymMax;
	char filename[128], line[1024], tempElem[8];
	FILE *input, *output;
	double temp, diff, diffnew, *C, Cdiff, Cdiffnew, C1, C2;
	VecR *Vdr1, *Vdr2, Vdr; 
	char **element1, **element2;

	//obtain representive clusters of goal big structure
	GoalClusterExtract ();

	//obtain the maximum symmetry value of representive cluster
	NsymMax = 0;
	for (i = 0; i < Nrc_goal; i ++) {
		if (rc_goal[i].Nsym > NsymMax) NsymMax = rc_goal[i].Nsym;
	}
	printf ("the maximum symmetry value of goal representive cluster is %d\n", NsymMax);

	//read initial atom structure of map small structure
	MapStructure ();

	//match the whole elements concentration between the big and small structure
	if (strcmp (cell_order, "disorder") == 0) {
		AllocMem (C, Nelems, double);
		for (i = 0; i < Nelems; i ++) {
			C[i] = 0;
			for (j = 0;j < Ngoal; j ++) {
				if (strcmp (goal[j].elem, elemName[i]) == 0) C[i] ++;
			}
			C[i] /= (double)Ngoal;
		}
		printf ("\nbig structure:");
		for (i = 0; i < Nelems; i ++) printf ("%s %.3f ", elemName[i], C[i]);
		printf ("\n");
		k = 0;
		for (i = 0; i < Nelems; i ++) {
			for (j = k; j < nMol; j ++) {
				sprintf (mol[j].elem, "%s", elemName[i]);
			}
			k += (int)(C[i] * (double)nMol + 0.5);
		}
		for (i = 0; i < Nelems; i ++) {
			C[i] = 0;
			for (j = 0;j < nMol; j ++) {
				if (strcmp (mol[j].elem, elemName[i]) == 0) C[i] ++;
			}
			C[i] /= (double)nMol;
		}
		printf ("small structure:");
		for (i = 0; i < Nelems; i ++) printf ("%s %.3f ", elemName[i], C[i]);
		printf ("\n\n");
		for (i = 0; i < nMol * nMol; i ++) { //random exchange
			while (1) {
				j = rand() % nMol;
				k = rand() % nMol;
				if (strcmp (mol[j].elem, mol[k].elem) != 0) break;
			}
			sprintf (tempElem, "%s", mol[j].elem);
			sprintf (mol[j].elem, "%s", mol[k].elem);
			sprintf (mol[k].elem, "%s", tempElem);
		}
	}
	PrintMapStructure ();
	system ("cp ../out/mapping.xyz ../in/trial.xyz");

	//build cluster cell of small structure
	for (i = 0; i < nMol; i ++) mol[i].NCEatom = 1;
	for (i = 0; i < nMol; i ++) {
		AllocMem (mol[i].CEatom, 1, int);
		mol[i].CEatom[0] = i;
		for (j = 0; j < nMol; j ++) {
			if (i == j) continue;
			VSub (Vdr, mol[j].r, mol[i].r);
			PeriodBoundary(Vdr);
			if (Vdr.x > -0.01 && Vdr.x  < regionCE.x - 0.01 &&
			    Vdr.y > -0.01 && Vdr.y  < regionCE.y - 0.01 &&
			    Vdr.z > -0.01 && Vdr.z  < regionCE.z - 0.01) {
				mol[i].NCEatom ++;
				ReAllocMem (mol[i].CEatom, mol[i].NCEatom, int);
				k = mol[i].NCEatom - 1;
				mol[i].CEatom[k] = j;
			}
		}
	}

	//allocate memory
	for (i = 0; i < nMol; i ++) {
		AllocMem (mol[i].Nelem, Nelems, int);
		AllocMem (mol[i].Celem, Nelems, double);
	}
	AllocMem (C_map, NC, CTABLE);  //the number of clusters for each atomic concentration
	for (i = 0; i < NC; i ++) AllocMem (C_map[i].Nc, Nelems, int);
	AllocMem (rc_map, nMol, Mol);
	for (i = 0; i < nMol; i ++) {
		AllocMem (rc_map[i].Nelem, Nelems, int);
		AllocMem (rc_map[i].Celem, Nelems, double);
		AllocMem (rc_map[i].CEatom, NC-1, int);
	}
	AllocMem (mat, (nMol+Ngoal), MATCH);

	loop = 0;
	diff = 1.e100;
	Cdiff = 1.e100;
	output = fopen ("../out/diff.dat", "w");
	fprintf (output, "task        Cdiff       CEdiff\n");
	printf ("task        Cdiff       CEdiff\n");
	while (1) {
		if (loop >= NLOOP) break;

		//random exchange the element between atoms A and B
		while (1) {
			A = rand() % nMol;
			B = rand() % nMol;
			if (strcmp (mol[A].elem, mol[B].elem) != 0) break;
		}
		sprintf (tempElem, "%s", mol[A].elem);
		sprintf (mol[A].elem, "%s", mol[B].elem);
		sprintf (mol[B].elem, "%s", tempElem);

		//obtain concentration table of small structure
		MapConcentrationTable ();
		Cdiffnew = 0;
		for (i = 0; i < NC; i ++) {
			for (j = 0; j < Nelems; j ++) {
				C1 = (double)(C_map[i].Nc[j]) / (double)nMol;
				C2 = (double)(C_goal[i].Nc[j]) / (double)Ngoal;
				if (error_map == 1) Cdiffnew += Sqr(C1 - C2);
				else if (error_map == 2) Cdiffnew += fabs(C1 - C2);
			}
		}

		if (loop < NLOOP / 2 && strcmp (cell_order, "disorder") == 0) {
			if (Cdiffnew > Cdiff && (double)rand() / RAND_MAX > 0.00001){
				sprintf (tempElem, "%s", mol[A].elem);
				sprintf (mol[A].elem, "%s", mol[B].elem);
				sprintf (mol[B].elem, "%s", tempElem);
			} else {
				Cdiff = Cdiffnew;
			}

			if (loop % (int)(NLOOP / 100) == 0) {
				fprintf (output, "%.3f%% %12.3e\n", (double)loop / (double)NLOOP * 100., Cdiff);
				printf ("%.3f%% %12.3e\n", (double)loop / (double)NLOOP * 100., Cdiff);

				PrintConcentration ();
				PrintMapStructure ();
			}
		} else {
			//obtain concentration table and representive clusters of small structure
			MapClusterExtract ();	

			//match the representive clusters between the big and small structure, then calculate the diff
			AllocMem (Vdr1, NC-1, VecR);
			AllocMem (Vdr2, NC-1, VecR);
			AllocMem2 (element1, NC-1, 8, char);
			AllocMem2 (element2, NC-1, 8, char);
			for (i = 0; i < Nrc_goal; i ++) rc_goal[i].match_flag = 0;
			for (i = 0; i < Nrc_map; i ++) rc_map[i].match_flag = 0;
			diffnew = 0.;
			Nmat = 0;
			for (i = 0; i < Nrc_map; i ++) {
				for (j = 0; j < Nrc_goal; j ++) {
					if (rc_goal[j].Nsym < NsymMax / 1000) continue;
					flag_elem = 1;
					for (k = 0; k < Nelems; k ++) {
						if (rc_map[i].Nelem[k] != rc_goal[j].Nelem[k]) flag_elem = 0;
					}

					if (flag_elem == 1) {
						//mol -> Vdr
						for (i1 = 0; i1 < rc_map[i].NCEatom; i1 ++) {
							k1 = rc_map[i].CEatom[i1];
							sprintf(element1[i1], "%s", mol[k1].elem);
							VSub (Vdr1[i1], mol[k1].r, rc_map[i].r);
							PeriodBoundary (Vdr1[i1]);
							PeriodBoundaryCE (Vdr1[i1]);
						}
						for (i2 = 0; i2 < rc_goal[j].NCEatom; i2 ++) {
							k2 = rc_goal[j].CEatom[i2];
							sprintf(element2[i2], "%s", goal[k2].elem);
							VSub (Vdr2[i2], goal[k2].r, rc_goal[j].r);
							PeriodBoundaryGoal (Vdr2[i2]);
							PeriodBoundaryCE (Vdr2[i2]);
						}
						flag_sym = Symmetry (Vdr1, element1, Vdr2, element2);
						if (flag_sym == 1) {
							rc_map[i].match_flag = 1;
							rc_goal[j].match_flag = 1;
							Nmat ++;
							mat[Nmat-1].small_rc_id = i;
							mat[Nmat-1].Csmall_rc = (double)(rc_map[i].Nsym) / (double)nMol;
							mat[Nmat-1].big_rc_id = j;
							mat[Nmat-1].Cbig_rc = (double)(rc_goal[j].Nsym) / (double)Ngoal;
							if (error_map == 1) {
								diffnew += Sqr((double)(rc_map[i].Nsym) / (double)nMol - (double)(rc_goal[j].Nsym) / (double)Ngoal);
							} else if (error_map == 2) {
								diffnew += fabs((double)(rc_map[i].Nsym) / (double)nMol - (double)(rc_goal[j].Nsym) / (double)Ngoal);
							}
						}
					}
				}
			}
			for (i = 0; i < Nrc_map; i ++) {
				if (rc_map[i].match_flag == 0) {
					Nmat ++;
					mat[Nmat-1].small_rc_id = i;
					mat[Nmat-1].Csmall_rc = (double)(rc_map[i].Nsym) / (double)nMol;
					mat[Nmat-1].big_rc_id = -1;
					mat[Nmat-1].Cbig_rc = 0.;
					if (error_map == 1) diffnew += Sqr((double)(rc_map[i].Nsym) / (double)nMol);
					else if (error_map == 2) diffnew += fabs((double)(rc_map[i].Nsym) / (double)nMol);
				}
			}
			for (i = 0; i < Nrc_goal; i ++) {
				if (rc_goal[i].match_flag == 0) {
					Nmat ++;
					mat[Nmat-1].small_rc_id = -1;
					mat[Nmat-1].Csmall_rc = 0.;
					mat[Nmat-1].big_rc_id = i;
					mat[Nmat-1].Cbig_rc = (double)(rc_goal[i].Nsym) / (double)Ngoal;
					if (error_map == 1) diffnew += Sqr ((double)(rc_goal[i].Nsym) / (double)Ngoal);
					else if (error_map == 2) diffnew += fabs ((double)(rc_goal[i].Nsym) / (double)Ngoal);
				}
			}
			if (diffnew > diff && (double)rand() / RAND_MAX > 0.00001){
				sprintf (tempElem, "%s", mol[A].elem);
				sprintf (mol[A].elem, "%s", mol[B].elem);
				sprintf (mol[B].elem, "%s", tempElem);
			} else {
				diff = diffnew;
				Cdiff = Cdiffnew;
			}

			//free memory
			free (Vdr1);
			free (Vdr2);
			free (element1[0]);
			free (element1);
			free (element2[0]);
			free (element2);

			if (loop % 10 == 0) {
				printf ("%.3f%% %12.3e %12.3e\n", (double)loop / (double)NLOOP * 100., Cdiff, diff);
				fprintf (output, "%.3f%% %12.3e %12.3e\n", (double)loop / (double)NLOOP * 100., Cdiff, diff);

				PrintRepresentiveCluster ();
				PrintMatch ();
				PrintConcentration ();
				PrintMapStructure ();
			}
		}

		loop ++;
	}
	fclose (output);

	//free memory
	for (i = 0; i < nMol; i ++) {
		free (mol[i].CEatom);
		free (mol[i].Nelem);
		free (mol[i].Celem);
	}
	for (i = 0; i < NC; i ++) free (C_map[i].Nc);
	free (C_map);
	for (i = 0; i < nMol; i ++) {
		free (rc_map[i].Nelem);
		free (rc_map[i].Celem);
		free (rc_map[i].CEatom);
	}
	free (rc_map);

	//free memory
	for (i = 0; i < Ngoal; i ++) {
		free (goal[i].Nelem);
		free (goal[i].Celem);
		free (goal[i].CEatom);
	}
	free (goal);
	for (i = 0; i < NC; i ++) free (C_goal[i].Nc);
	free (C_goal);
	for (i = 0; i < Ngoal; i ++) {
		free (rc_goal[i].Nelem);
		free (rc_goal[i].Celem);
		free (rc_goal[i].CEatom);
	}
	free (rc_goal);

	free (C);

	free (mat);
}

void MapStructure ()
{
	int i;
	FILE *input;
	char line[1024];
	double temp;

	input = fopen ("../in/trial.xyz", "r");
	fgets (line, 1024, input);
	sscanf (line, "%d", &nMol);
	fgets (line, 1024, input);
	sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", &region.x, &temp, &temp, &temp, &region.y, &temp, &temp, &temp, &region.z);
	printf ("the small structure region:[%f %f %f]\n", region.x, region.y, region.z); 
	AllocMem (mol, nMol, Mol);
	for (i = 0; i < nMol; i ++) {
		fgets (line, 1024, input);
		sscanf (line, "%s %lg %lg %lg", mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z);
	}
}

void MapConcentrationTable ()
{
	int i, j, k, id;
	VecR Vdr;
	FILE *output;
	char filename[128];

	//calculate atomic concentration of each cluster
	for (i = 0; i < nMol; i ++) {
		for (j = 0; j < Nelems; j ++) mol[i].Nelem[j] = 0.;
		for (k = 0; k < mol[i].NCEatom; k ++) {
			id = mol[i].CEatom[k];
			for (j = 0; j < Nelems; j ++) {
				if (strcmp(mol[id].elem, elemName[j]) == 0) mol[i].Nelem[j] ++;
			}
		}
		for (j = 0; j < Nelems; j ++) mol[i].Celem[j] = (double)(mol[i].Nelem[j]) / (double)(mol[i].NCEatom);
	}

	//obtain concentration table of small structure
	for (i = 0; i < nMol; i ++) {
		if (mol[i].NCEatom + 1 != NC) printf ("warning: box region may be wrong\n");
	}

	//set initial value
	for (i = 0; i < NC; i ++) {
		C_map[i].C = (double)i / (double)(NC - 1);
		for (j = 0; j < Nelems; j ++) C_map[i].Nc[j] = 0;
	}

	//obtain concentration table
	for (i = 0; i < nMol; i ++) {
		for (j = 0; j < Nelems; j ++) {
			k = mol[i].Nelem[j];
			C_map[k].Nc[j] ++;
		}
	}
}

void MapClusterExtract ()
{
	MapConcentrationTable ();
	Nrc_map = ClusterExtract (nMol, mol, C_map, rc_map, 2); //obtain concentration table and representive clusters of small structure
//	printf ("the representive cluster number: %d\n", Nrc_map);
}

void GoalClusterExtract ()
{
	FILE *input, *output;
	char line[1024], filename[128], command[128];
	int N, i, j, k, l, flag = 0, id;
	double temp;

	VecR Vdr;

	//read the region of extrected clusters cell
	sprintf (filename, "../in/CONTCAR");
	if ((input = fopen (filename, "r")) == 0) {
		printf ("%s doesn't exist\n", filename);
		exit (1);
	}
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &regionCE.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &regionCE.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &regionCE.z);
	VMul (regionCE, regionCE, clusterSize);
	printf ("the cluster region: [%.3f %.3f %.3f]\n", regionCE.x, regionCE.y, regionCE.z);
	fclose (input);

	//read goal atom structure
	sprintf (filename, "../in/benchmark.xyz");
	printf ("reading atom information from benchmark.xyz ...\n");
	input = fopen (filename, "r");

	fgets (line, 1024, input);
	sscanf (line, "%d", &Ngoal);
	fgets (line, 1024, input);
	sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", &region_goal.x, &temp, &temp, &temp, &region_goal.y, &temp, &temp, &temp, &region_goal.z);
	printf ("the benchmark supercell region:[%f %f %f]\n", region_goal.x, region_goal.y, region_goal.z); 
	printf ("building neighbour list...\n");
	AllocMem (goal, Ngoal, Mol);
	for (i = 0; i < Ngoal; i ++) {
		fgets (line, 1024, input);
		sscanf (line, "%s %lg %lg %lg", goal[i].elem, &goal[i].r.x, &goal[i].r.y, &goal[i].r.z);
		if (strcmp(goal[i].elem, "vacancy") == 0) sprintf (goal[i].elem, "%s", elemName[0]);
	}

	//build cluster cell
	for (i = 0; i < Ngoal; i ++) goal[i].NCEatom = 1;
	#pragma omp parallel for private(j, Vdr, k) num_threads(nthreads)
	for (i = 0; i < Ngoal; i ++) {
		AllocMem (goal[i].CEatom, 1, int);
		goal[i].CEatom[0] = i;
		for (j = 0; j < Ngoal; j ++) {
			if (i == j) continue;
			VSub (Vdr, goal[j].r, goal[i].r);
			PeriodBoundaryGoal(Vdr);
			if (Vdr.x > -0.01 && Vdr.x  < regionCE.x - 0.01 &&
			    Vdr.y > -0.01 && Vdr.y  < regionCE.y - 0.01 &&
			    Vdr.z > -0.01 && Vdr.z  < regionCE.z - 0.01) {
				goal[i].NCEatom ++;
				ReAllocMem (goal[i].CEatom, goal[i].NCEatom, int);
				k = goal[i].NCEatom - 1;
				goal[i].CEatom[k] = j;
			}
		}
	}

	//calculate atomic concentration of each cluster
	for (i = 0; i < Ngoal; i ++) {
		AllocMem (goal[i].Nelem, Nelems, int);
		AllocMem (goal[i].Celem, Nelems, double);
	}
	for (i = 0; i < Ngoal; i ++) {
		for (j = 0; j < Nelems; j ++) goal[i].Nelem[j] = 0.;
		for (k = 0; k < goal[i].NCEatom; k ++) {
			id = goal[i].CEatom[k];
			for (j = 0; j < Nelems; j ++) {
				if (strcmp(goal[id].elem, elemName[j]) == 0) goal[i].Nelem[j] ++;
			}
		}
		for (j = 0; j < Nelems; j ++) goal[i].Celem[j] = (double)(goal[i].Nelem[j]) / (double)(goal[i].NCEatom);
	}

	NC = goal[0].NCEatom + 1; //the number of atomic concentration
	for (i = 1; i < Ngoal; i ++) {
		if (goal[i].NCEatom + 1 != NC) printf ("warning: box region may be wrong\n");
	}
	AllocMem (C_goal, NC, CTABLE);  //the number of clusters for each atomic concentration
	for (i = 0; i < NC; i ++) AllocMem (C_goal[i].Nc, Nelems, int);
	AllocMem (rc_goal, Ngoal, Mol);
	for (i = 0; i < Ngoal; i ++) {
		AllocMem (rc_goal[i].Nelem, Nelems, int);
		AllocMem (rc_goal[i].Celem, Nelems, double);
		AllocMem (rc_goal[i].CEatom, goal[i].NCEatom, int);
	}

	//set initial value
	for (i = 0; i < NC; i ++) {
		C_goal[i].C = (double)i / (double)(NC - 1);
		for (j = 0; j < Nelems; j ++) C_goal[i].Nc[j] = 0;
	}

	//obtain concentration table
	for (i = 0; i < Ngoal; i ++) {
		for (j = 0; j < Nelems; j ++) {
			k = goal[i].Nelem[j];
			C_goal[k].Nc[j] ++;
		}
	}
//	for (i = 0; i < NC; i ++) {
//		for (j = 0; j < Nelems; j ++) printf ("%.3f (%s %d) ", C[i].C, elemName[j], C[i].Nc[j]);
//		printf ("\n");
//	}

	//obtain representive clusters of goal supercell
	Nrc_goal = ClusterExtract (Ngoal, goal, C_goal, rc_goal, 1);
	printf ("the representive cluster number: %d\n", Nrc_goal);

	//print concentration
	output = fopen ("../out/CE_benchmark/concentration.dat", "w");
	printf ("\nC          ");
	for (i = 0; i < Nelems; i ++) printf ("N(%s)     proportion(%s) ", elemName[i], elemName[i]);
	printf ("\n");
	fprintf (output, "C          ");
	for (i = 0; i < Nelems; i ++) fprintf (output, "N(%s)      proportion(%s) ", elemName[i], elemName[i]);
	fprintf (output, "\n");
	for (i = 0; i < NC; i ++) {
		printf ("%.3f ", C_goal[i].C);
		for (j = 0; j < Nelems; j ++) printf ("%10d %12.3f ", C_goal[i].Nc[j], (double)(C_goal[i].Nc[j]) / (double)Ngoal);
		printf ("\n");
		fprintf (output, "%.3f ", C_goal[i].C);
		for (j = 0; j < Nelems; j ++) fprintf (output, "%10d %12.3f ", C_goal[i].Nc[j], (double)(C_goal[i].Nc[j]) / (double)Ngoal);
		fprintf (output, "\n");
	}
	fclose (output);

	//print representive clusters
	for (i = 0; i < Nrc_goal; i ++) {
		sprintf (filename, "../out/CE_benchmark/CE%d.xyz", i);
		output = fopen (filename, "w");
		fprintf (output, "%d\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" Properties=species:S:1:pos:R:3\n", rc_goal[i].NCEatom, regionCE.x, regionCE.y, regionCE.z);
		for (j = 0; j < rc_goal[i].NCEatom; j ++) {
			id = rc_goal[i].CEatom[j];
			VSub (Vdr, goal[id].r, rc_goal[i].r);
			PeriodBoundaryGoal(Vdr);
			fprintf (output, "%s %e %e %e\n", goal[id].elem, Vdr.x, Vdr.y, Vdr.z);
		}
		fclose (output);
	}

	sprintf (filename, "../out/CE_benchmark/CE.dat");
	output = fopen (filename, "w");
	fprintf (output, "CE N\n");
	for (i = 0; i < Nrc_goal; i ++) fprintf (output, "%d %d\n", i, rc_goal[i].Nsym);
	fclose (output);
}

int ClusterExtract (int N, Mol *cluster, CTABLE *C, Mol *rc, int mode) //obtain concentration table and representive clusters
{
	int i, j, k, i1, j1, k1, i2, j2, k2, Nrc, flag_add, flag_elem, flag_sym;
	VecR *Vdr1, *Vdr2; 
	char **element1, **element2;

	//obtain representive clusters
	AllocMem (Vdr1, NC-1, VecR);
	AllocMem (Vdr2, NC-1, VecR);
	AllocMem2 (element1, NC-1, 8, char);
	AllocMem2 (element2, NC-1, 8, char);
	Nrc = 1;
	VCopy(rc[0].r, cluster[0].r); //set first representive clusters
	sprintf (rc[0].elem, "%s", cluster[0].elem);
	for (i = 0; i < Nelems; i ++) {
		rc[0].Nelem[i] = cluster[0].Nelem[i];
		rc[0].Celem[i] = cluster[0].Celem[i];
	}
	rc[0].NCEatom = cluster[0].NCEatom;
	for (i = 0; i < rc[0].NCEatom; i ++) {
		rc[0].CEatom[i] = cluster[0].CEatom[i];
	}
	rc[0].Nsym = 1;
	for (i = 1; i < N; i ++) {
		flag_add = 1; //if adding cluster to lib
		for (j = 0; j < Nrc; j ++) {
			flag_elem = 1;
			for (k = 0; k < Nelems; k ++) {
				if (cluster[i].Nelem[k] != rc[j].Nelem[k]) flag_elem = 0;
			}
			if (flag_elem == 1) {
				//mol -> Vdr
				for (i1 = 0; i1 < cluster[i].NCEatom; i1 ++) {
					k1 = cluster[i].CEatom[i1];
					if (mode == 1) {
						sprintf(element1[i1], "%s", goal[k1].elem);
						VSub (Vdr1[i1], goal[k1].r, cluster[i].r);
						PeriodBoundaryGoal (Vdr1[i1]);
					} else if (mode == 2) {
						sprintf(element1[i1], "%s", mol[k1].elem);
						VSub (Vdr1[i1], mol[k1].r, cluster[i].r);
						PeriodBoundary (Vdr1[i1]);
					}
					PeriodBoundaryCE (Vdr1[i1]);
				}
				for (i2 = 0; i2 < rc[j].NCEatom; i2 ++) {
					k2 = rc[j].CEatom[i2];
					if (mode == 1) {
						sprintf(element2[i2], "%s", goal[k2].elem);
						VSub (Vdr2[i2], goal[k2].r, rc[j].r);
						PeriodBoundaryGoal (Vdr2[i2]);
					} else if (mode == 2) {
						sprintf(element2[i2], "%s", mol[k2].elem);
						VSub (Vdr2[i2], mol[k2].r, rc[j].r);
						PeriodBoundary (Vdr2[i2]);
					}
					PeriodBoundaryCE (Vdr2[i2]);
				}

				flag_sym = Symmetry (Vdr1, element1, Vdr2, element2);

				if (flag_sym == 1) {
					flag_add = 0;
					rc[j].Nsym ++;
				}
			}
		}
		if (flag_add == 1) {
			Nrc ++;

			//add representive clusters
			VCopy(rc[Nrc-1].r, cluster[i].r);
			sprintf (rc[Nrc-1].elem, "%s", cluster[i].elem);
			for (j = 0; j < Nelems; j ++) {
				rc[Nrc-1].Nelem[j] = cluster[i].Nelem[j];
				rc[Nrc-1].Celem[j] = cluster[i].Celem[j];
			}
			rc[Nrc-1].NCEatom = cluster[i].NCEatom;
			for (j = 0; j < rc[Nrc-1].NCEatom; j ++) {
				rc[Nrc-1].CEatom[j] = cluster[i].CEatom[j];
			}
			rc[Nrc-1].Nsym = 1;
		}
		if (mode == 1 && i % (N / 100) == 0) printf ("representive clusters: %d, %.0f%%\n", Nrc, i * 100. / N);
	}

	free (Vdr1);
	free (Vdr2);
	free (element1[0]);
	free (element1);
	free (element2[0]);
	free (element2);

	return (Nrc);
}

int Symmetry (VecR *Vdr1, char **element1, VecR *Vdr2, char **element2)
{
	int flag, flag_find, i, j, k, l, k1, k2, Ndr;
	VecR Vdr, Vdrr, Vdr1s, Vdr2s, VdrNebr[1000], VdrMove;
	double dr, temp, x, y, z;
	FILE *output;
	char filename[128];

	//find the the move vectors
	VSub (Vdr, Vdr1[1], Vdr1[0]);
	PeriodBoundaryCEdr (Vdr);
	VCopy (VdrNebr[0], Vdr);
	Ndr = 1;
	for (i = 2; i < NC-1; i ++) {
		VSub (Vdr, Vdr1[i], Vdr1[0]);
		PeriodBoundaryCEdr (Vdr);
		flag = 0;
		for (j = 0; j < Ndr; j ++) {
			VSub (Vdrr, Vdr, VdrNebr[j]);
			dr = VLen (Vdrr);
			if (dr < 0.01) flag = 1;
		}
		if (flag == 0) {
			Ndr ++;
			VCopy (VdrNebr[Ndr-1], Vdr);
		}
	}
/*	
	printf ("move vectors: %d\n", Ndr);
	for (i = 0; i < Ndr; i ++) printf ("%.3f %.3f %.3f\n", VdrNebr[i].x, VdrNebr[i].y, VdrNebr[i].z);
	exit (1);
*/
	for (l = -1; l < Ndr; l ++) {
		for (k = 0; k < Nsym; k ++) {
			flag = 1;
			for (i = 0; i < NC-1; i ++) {
				//move
				if (l == -1) VCopy (VdrMove, Vdr1[i]);
				else VAdd (VdrMove, Vdr1[i], VdrNebr[l]);
				PeriodBoundaryCE (VdrMove);
				// symmetry
				x = VdrMove.x * sym[k][0][0] + VdrMove.y * sym[k][1][0] + VdrMove.z * sym[k][2][0];
				y = VdrMove.x * sym[k][0][1] + VdrMove.y * sym[k][1][1] + VdrMove.z * sym[k][2][1];
				z = VdrMove.x * sym[k][0][2] + VdrMove.y * sym[k][1][2] + VdrMove.z * sym[k][2][2];
				VSet (Vdr1s, x, y, z);
				PeriodBoundaryCE (Vdr1s);
				flag_find = 0;
				for (j = 0; j < NC-1; j ++) {
					VSub (Vdr, Vdr1s, Vdr2[j]);
					dr = VLen (Vdr);
					if (fabs(dr) < 0.01) {
						flag_find = 1;
						if (strcmp (element1[i], element2[j]) != 0) flag = 0;
					}
				}
				if (flag_find == 0) {
					printf ("atom: [%.3f %.3f %.3f]\n", Vdr1s.x, Vdr1s.y, Vdr1s.z);
					PrintCluster (NC-1, Vdr1, element1);
					PrintCluster (NC-1, Vdr2, element2);
					printf ("error: does not find the corresponding atoms when obtaining representive cluster\n");
					exit (1);
				}
			}
			if (flag == 1) return (1);
		}
	}

	return (0);
}

void PrintCluster (int N, VecR *c, char **element)
{
	FILE *output;
	char filename[128];
	int i;

	sprintf (filename, "../out/cluster.movie");
	output = fopen (filename, "a+");
	fprintf (output, "%d\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" Properties=species:S:1:pos:R:3\n", N, regionCE.x, regionCE.y, regionCE.z);
	for (i = 0; i < N; i ++) {
		fprintf (output, "%s %f %f %f\n", element[i], c[i].x, c[i].y, c[i].z);
	}

	fclose (output);
}

void PrintRepresentiveCluster ()
{
	int i, j, id;
	VecR Vdr;
	char filename[128];
	FILE *output;

	//print representive clusters
	system ("rm -rf ../out/CE_trial/CE*.xyz");
	for (i = 0; i < Nrc_map; i ++) {
		sprintf (filename, "../out/CE_trial/CE%d.xyz", i);
		output = fopen (filename, "w");
		fprintf (output, "%d\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" Properties=species:S:1:pos:R:3\n", rc_map[i].NCEatom, regionCE.x, regionCE.y, regionCE.z);
		for (j = 0; j < rc_map[i].NCEatom; j ++) {
			id = rc_map[i].CEatom[j];
			VSub (Vdr, mol[id].r, rc_map[i].r);
			PeriodBoundary(Vdr);
			fprintf (output, "%s %e %e %e\n", mol[id].elem, Vdr.x, Vdr.y, Vdr.z);
		}
		fclose (output);
	}
	sprintf (filename, "../out/CE_trial/CE.dat");
	output = fopen (filename, "w");
	fprintf (output, "CE N\n");
	for (i = 0; i < Nrc_map; i ++) fprintf (output, "%d %d\n", i, rc_map[i].Nsym);
	fclose (output);
}

void PrintConcentration ()
{
	int i, j;
	FILE *output;

	//print concentration
	output = fopen ("../out/CE_trial/concentration.dat", "w");
	fprintf (output, "C          ");
	for (i = 0; i < Nelems; i ++) fprintf (output, "N(%s)      proportion(%s) ", elemName[i], elemName[i]);
	fprintf (output, "\n");
	for (i = 0; i < NC; i ++) {
		fprintf (output, "%.3f ", C_map[i].C);
		for (j = 0; j < Nelems; j ++) fprintf (output, "%10d %12.3f ", C_map[i].Nc[j], (double)(C_map[i].Nc[j]) / (double)nMol);
		fprintf (output, "\n");
	}
	fclose (output);
}

void PrintMapStructure ()
{
	FILE *output;
	int n;
	
	output = fopen ("../out/mapping.xyz", "w");
	fprintf (output, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=species:S:1:pos:R:3\n", \
		 nMol, region.x, region.y, region.z);
	DO_MOL fprintf (output, "%s %f %f %f\n", mol[n].elem, mol[n].r.x, mol[n].r.y, mol[n].r.z); //--A
	fclose (output);
}

void PrintMatch ()
{
	FILE *output;
	int i;
	double diff;

	output = fopen ("../out/mapping.dat", "w");
	fprintf (output, "trial_id     C     benchmark_id    C         diff\n");
	for (i = 0; i < Nmat; i ++) {
		if (error_map == 1) diff = Sqr(mat[i].Csmall_rc - mat[i].Cbig_rc);
		else if (error_map == 2) diff = fabs(mat[i].Csmall_rc - mat[i].Cbig_rc);
		
		if (mat[i].small_rc_id == -1) {
			fprintf (output, " /    %10.3e %10d %10.3e %10.3e\n", mat[i].Csmall_rc, mat[i].big_rc_id, mat[i].Cbig_rc, diff);
		} else if (mat[i].big_rc_id == -1) {
			fprintf (output, "%3d   %10.3e         /  %10.3e %10.3e\n", mat[i].small_rc_id, mat[i].Csmall_rc, mat[i].Cbig_rc, diff);
		} else {
			fprintf (output, "%3d   %10.3e %10d %10.3e %10.3e\n", mat[i].small_rc_id, mat[i].Csmall_rc, mat[i].big_rc_id, mat[i].Cbig_rc, diff);
		}
	}

	fclose (output);
}
