#include "measurement.h"

int flag_SRO = 0;

//RDF
double **histRdf, rangeRdf;
int countRdf, limitRdf, sizeHistRdf, stepRdf;
double *ith_nebrR;
double *C_threshold;
char **elem_preci;

// alphai = 1 - PiAB / cB, cB is the concentration of atoms of B sort.
// PiAB = NiAB / NiA, PiAB is the probability to find a B atom on the i-th sphere of coordination of an A atom located at the origin of coordinates.
double ShortRangeOrder (char *A, char *B, char *filename)
{
	int n, m, nebr, i;
	int NiAB, NiA, NA;
	double alphai, piAB, cB, PiAB, ci;
	FILE *output;

	NiAB = 0;
	NiA = 0;
	NA = 0;
	i = 0;
	#pragma omp parallel for reduction (+:NA, NiA, NiAB) private (m, nebr) num_threads (nthreads)
	DO_ATOM {
		if (strcmp(atom[n].element, A) == 0) {
			NA ++;
			NiA += atom[n].nebrLists[i].Nnebrs;
			for (m = 0; m < atom[n].nebrLists[i].Nnebrs; m++) {
				nebr = atom[n].nebrLists[i].nebr[m];
				if (strcmp(atom[nebr].element, B) == 0) NiAB ++;
			}
		} 
	}
	PiAB = (double)NiAB / (double)NiA;
	cB = 1. - (double)NA / (double)Natoms[0];
	alphai = 1. - PiAB / cB;

	printf ("SRO: %f, ", alphai);

	output = WriteFile (filename);
	if (flag_SRO == 0) fprintf (output, "step SRO\n");
	fprintf (output, "%lld %f\n", stepCount, alphai);

	flag_SRO = 1;
	fclose (output);

	return (alphai);
}

void MeasurementRdf ()
{
	double rb, dx, *integral;
	int n, i, j;
	char filename[128];
	FILE *output, *output1;

	sprintf (filename, "../out/RDF.dat");
	output = WriteFile (filename);
	sprintf (filename, "../out/RDF_integral.dat");
	output1 = WriteFile (filename);
	fprintf (output, "dr(A) rdf ");
	fprintf (output1, "dr(A) rdf ");
	for (i = 0; i < Nelems; i ++) {
		for (j = 0; j < Nelems; j ++) {
			fprintf (output, "%s(%s) ", elemName[i], elemName[j]);
			fprintf (output1, "%s(%s) ", elemName[i], elemName[j]);
		}
	}
	fprintf (output, "\n");
	fprintf (output1, "\n");
	
	dx = rangeRdf / sizeHistRdf;
	AllocMem (integral, Nelems*Nelems+1, double);
	for (i = 0; i < Nelems*Nelems+1; i ++) integral[i] = 0.;
	for (n = 0; n < sizeHistRdf; n ++){
		rb = (n + 0.5) * rangeRdf / sizeHistRdf; //--dimensionless
		fprintf (output, "%f ", rb*lUnit*1.e10);
		fprintf (output1, "%f ", rb*lUnit*1.e10);
		for (i = 0; i < Nelems*Nelems+1; i ++) {
			fprintf (output, "%f ", histRdf[i][n]);
			integral[i] += histRdf[i][n] * dx;
			fprintf (output1, "%f ", integral[i]);
		}
		fprintf (output, "\n");
		fprintf (output1, "\n");
	}

/*
	//extremum
	int i, testnu = sizeHistRdf / 10 + 2, ifmax = 1, ifmin = 1;
	for (n = testnu; n < sizeHistRdf - testnu; n++){
		for (i = 1; i <= testnu; i ++){
			if (histRdf[n] <= histRdf[n-i] || histRdf[n] <= histRdf[n+i]) ifmax = 0;
			if (histRdf[n] >= histRdf[n-i] || histRdf[n] >= histRdf[n+i]) ifmin = 0;
		}
		if (ifmax == 1){ 
//			printf("maxr = %.4f, maxrdf = %.4f, ", (n + 0.5) * dx, histRdf[n]);
		}
		if (ifmin == 1){
//			printf("minr = %.4f, minrdf = %.4f, ", (n + 0.5) * dx, histRdf[n]);
			rRdf = (n + 0.5) * dx;
			break;
		}
		ifmax = 1;
		ifmin = 1;
	}
//	printf("\n\n");
*/
	free (integral);
	fclose (output);
	fclose (output1);
}

void EvalRdf ()
{
	VecR dr;
	double deltaR, *normFac, rr;
	int j1, j2, n, k, nthreads, NRdf, i, j, count;

	NRdf = Nelems*Nelems+1.;
	if (countRdf == 0) {
		for (k = 0; k < NRdf; k ++)
			for (n = 0; n < sizeHistRdf; n ++) histRdf[k][n] = 0.;
	}
	if (rangeRdf > rCut + rNebrShell) {
		printf ("error(measurement.c): rangeRdf is too large\n");
		exit (1);
	}
	deltaR = rangeRdf / sizeHistRdf;

	nthreads = 22;
	#pragma omp parallel for private (j1,j2,dr,rr,n,i,j) num_threads (nthreads)
	for (k = 0; k < nebrTabLen; k ++) {
		j1 = nebrTab[2 * k];
		j2 = nebrTab[2 * k + 1];
		VSub (dr, mol[j1].r, mol[j2].r);
		VWrapAll (dr);
		rr = VLenSq (dr);
		if (rr < Sqr (rangeRdf)) {
			n = sqrt (rr) / deltaR;
			#pragma omp atomic
			++ histRdf[0][n];
			for (i = 0; i < Nelems; i ++) {
				for (j = 0; j < Nelems; j ++) {
					if ((strcmp(elemName[i], mol[j1].elem) == 0 && \
					     strcmp(elemName[j], mol[j2].elem) == 0) || \
					    (strcmp(elemName[i], mol[j2].elem) == 0 && \
					     strcmp(elemName[j], mol[j1].elem) == 0)) {
						#pragma omp atomic
						++ histRdf[j+i*Nelems+1][n];
					}
				}
			}
		}
	}

	++ countRdf;
	if (countRdf == limitRdf) {
		AllocMem (normFac, Nelems+1, double);
		for (n = 0; n <= Nelems; n ++)
			normFac[n] = VProd (region) / (2. * M_PI * Cube (deltaR) * Sqr (Natoms[n]) * countRdf);
		for (n = 0; n < sizeHistRdf; n ++) {
			histRdf[0][n] *= normFac[0] / Sqr (n - 0.5);
			for (i = 0; i < Nelems; i ++) {
				for (j = 0; j < Nelems; j ++) {
					histRdf[j+i*Nelems+1][n] *= normFac[i+1] / Sqr (n - 0.5);
				}
			}
		}
		free (normFac);
		MeasurementRdf ();
		countRdf = 0;
	}
}

void IdenPrecipitates ()
{
	FILE *input, *output;
	char filename[128], line[1024], title[1024];
	int m, n, N, id, temp, i, j, k, l, j1, j2, Nmovie, N_preci, Ncore, flag, nebr_flag = 0, preci_flag = 0;
	double r, rMin, tmp, V, V0, rnebr, d, densi, C, ith_nebrR_max, *drmin;
	VecR dr, *coreR, Vdr;
	int element_id, Nnebr, Cmax_id, *core_flag, preci_id, *Natoms_preci, d_profile[100];

	//Unit
	TUnit = 300.; //--K
	EUnit = TUnit * kB; //--J
	VUnit = EUnit / eleChar; //--eV
	lUnit = 1.e-10; //--m
	eUnit = EUnit; //--J

	sprintf (filename, "../out/md/out.movie");
	input = ReadFile (filename);

	//read out.movie
	i = nMol = Nmovie = 0;
	system ("rm -f ../out/precipitates*.movie ../out/precipitates*.dat");
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;

		// read data
		if (i == nMol) {
			Nmovie ++;
			sscanf (line, "%d", &nMol);
			if (nebr_flag == 0) {
				AllocMem (mol, nMol, Mol);
				AllocMem (coreR, nMol, VecR);
				AllocMem (core_flag, nMol, int);
			}
			fgets (line, 1024, input);
			sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", \
				&region.x, &tmp, &tmp, &tmp, &region.y, &tmp, &tmp, &tmp, &region.z); //--A
			if (isnormal(VProd (region)) == 0) {
				printf ("error(measurement.c): out.movie format, region(%f %f %f)\n", region.x, region.y, region.z);
				exit (1);
			}
			VScale (region, 1. / (lUnit * 1.e10)); //dimensionless
			fgets (line, 1024, input);		
			i = 0;
		}
		sscanf (line, "%d %s %lg %lg %lg", &temp, mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z); //--A
		VScale (mol[i].r, 1. / (lUnit * 1.e10)); //dimensionless

		/*******************identify precipitates*********************/
		if (i == nMol - 1) {
			if (nebr_flag == 0) {
				//set initial parameter and allocate memory, build neighbour list
				rCut = 10.; // --A
				rCut /= lUnit * 1.e10; //dimensionless
				nebrTabMax = nebrTabFac * nMol;
				AllocMem (nebrTab, 2 * nebrTabMax, int);
				DO_MOL AllocMem (mol[n].nebr_id, nebrTabFac, int);
				printf ("\n-- movie: %d, build neighbour list...\n", Nmovie);
				VSCopy (cells, 1. / (rCut + rNebrShell), region);
				AllocMem (cellList, VProd (cells) + nMol, int);
				BuildNebrList ();
				free (cellList);

				//find ith_nebrR
				ith_nebrR_max = 0;
				for (n = 0; n < Nelems; n ++) ith_nebrR_max = Max (ith_nebrR_max, ith_nebrR[n]);
				AllocMem (drmin, (int)ith_nebrR_max+1, double);
				drmin[0] = 0.;
				for (m = 1; m < (int)ith_nebrR_max+1; m ++) drmin[m] = 1.e10;
				for (m = 1; m < (int)ith_nebrR_max+1; m ++) {
					for (n = 1; n < nMol; n ++) {
						j1 = 0;
						j2 = n;
						VSub (Vdr, mol[j1].r, mol[j2].r);
						PeriodBoundary (Vdr);
						r = VLen (Vdr);
						if ((drmin[m]) > r && r > drmin[m-1]+0.01) drmin[m] = r;	
					}
					printf ("the %d nearest proximity distance is %f\n", m, drmin[m]);
				}
				rnebr = drmin[(int)ith_nebrR_max];
				printf ("the cut radius is %f\n", rnebr);
				free (drmin);
				VSCopy (cells, 1. / (rnebr), region);
				AllocMem (cellList, VProd (cells) + nMol, int);
				GetNebrs (rnebr);
				free (cellList);
				nebr_flag = 1;
			}
			
			printf ("\n-- movie: %d, identify precipitates...\n", Nmovie);
			DO_MOL mol[n].interface_flag = 1;
			for (element_id = 0; element_id < 1; element_id ++) {
				DO_MOL {
					N = 0;
					Nnebr = 0;
					for (k = 0; k < mol[n].Nnebr; k ++) {
						id = mol[n].nebr_id[k];
						VSub (dr, mol[n].r, mol[id].r);
						VWrapAll (dr);
						r = VLen (dr);
						if (strcmp(mol[id].elem, elem_preci[element_id]) == 0 && \
						    r < rMin * ith_nebrR[element_id]) N ++;
						if (r < rMin * ith_nebrR[element_id]) Nnebr ++;
					}
					if (strcmp(mol[n].elem, elem_preci[element_id]) == 0) N ++;
					mol[n].C = (double)N / (double)(Nnebr + 1);
				}

				DO_MOL {
					mol[n].preci_flag = 0;
					VZero (coreR[n]);
				}
				N_preci = 0; //the number of precipetitions
				Ncore = 0; //the number of the cores of precipitations shells
				DO_MOL {
					if (mol[n].C > C_threshold[element_id]) {
						//identify the id of precipetitions (core_flag)
						VCopy (coreR[Ncore], mol[n].r);
						flag = 1;
						if (Ncore > 0) {
							preci_id = -1;
							for (j = 0; j < Ncore; j ++) {
								VSub (dr, coreR[Ncore], coreR[j]);
								VWrapAll (dr);
								r = VLen (dr);
								if (r <  2. * rMin * ith_nebrR[element_id]) {
									flag = 0;
									if (preci_id == -1) {
										preci_id = core_flag[j];
										core_flag[Ncore] = preci_id;
									} else {
										if (core_flag[j] != preci_id) {
											for (k = 0; k < Ncore; k ++) {
												if (core_flag[k] == core_flag[j]) core_flag[k] = preci_id;
											}
										}
									}
								}
							}
						}
						if (flag == 1) {
							N_preci ++;
							core_flag[Ncore] = N_preci;
						}
						Ncore ++;
						//identify the precitations atoms
						mol[n].preci_flag = 1;
						for (k = 0; k < mol[n].Nnebr; k ++) {
							id = mol[n].nebr_id[k];
							VSub (dr, mol[n].r, mol[id].r);
							VWrapAll (dr);
							r = VLen (dr);
							if (r < rMin * ith_nebrR[element_id]) {
								mol[id].preci_flag = 1;
							}
						}
					}
				}
				//identify the number of precipetitions
				k = 1;
				l = 0;
				for (k = 1; k <= N_preci; k ++) {
					flag = 0;
					for (j = 0; j < Ncore; j ++) {
						if (core_flag[j] == k) flag = 1;
					}
					if (flag == 1) l ++;
				}
				N_preci = l;

				//the number of precitations atoms
				N = 0;
				DO_MOL {
					if (mol[n].preci_flag == 1) N ++;
				}

				//mark the precitations atoms belong to which core
				for (j = 0; j < Ncore; j ++) {
					DO_MOL {
						if (mol[n].preci_flag != 1) continue;
						VSub (dr, mol[n].r, coreR[j]);
						VWrapAll (dr);
						r = VLen (dr);
						if (r < rMin * ith_nebrR[element_id]) mol[n].id = core_flag[j];
					}
				}
				//obtain the size of each precipetition
				AllocMem (Natoms_preci, Ncore, int);
				for (j = 0; j < Ncore; j ++) Natoms_preci[j] = 0;
				DO_MOL {
					if (mol[n].preci_flag == 1) Natoms_preci[mol[n].id] ++;
				}
				for (k = 0; k < 100; k ++) d_profile[k] = 0;
				for (j = 0; j < Ncore; j ++) {
					if (Natoms_preci[j] != 0) {
						V0 = VProd (region) / nMol; //dimensionless
						d = pow (6. * (V0 * Natoms_preci[j]) / M_PI, 1. / 3.); //dimensionless
						d *= lUnit * 1.e9; //--nm
						for (k = 0; k < 100; k ++) {
							if (d > (double)k * 0.1 && d < (double)(k + 1) * 0.1) d_profile[k] ++;
						}
					}
				}
				sprintf (filename, "../out/d_profile_%s.dat", elem_preci[element_id]);
				output = WriteFile (filename);
				fprintf (output, "d(nm) density(m-3) PSD\n");
				for (k = 0; k < 100; k ++) fprintf (output, "%f %e %d\n", (double)k * 0.1+0.05, (double)d_profile[k] / (VProd(region) * Cube(lUnit)), d_profile[k]);
				free (Natoms_preci);
				fclose (output);

/*
				//identify the number of precipetitions
				N_preci = 0;
				DO_MOL {
					if (mol[n].preci_flag == 1) {
						Cmax_id = -1;
						for (k = 0; k < mol[n].Nnebr; k ++) {
							id = mol[n].nebr_id[k];
							VSub (dr, mol[n].r, mol[id].r);
							VWrapAll (dr);
							r = VLen (dr);
							if (mol[n].C < mol[id].C && r < rMin * ith_nebrR[element_id]) Cmax_id = id;
						}
						if (Cmax_id == -1) N_preci ++;
					}
				}
*/

				//output precipitation file
				sprintf (filename, "../out/precipitates_%s.dat", elem_preci[element_id]);
				output = WriteFile (filename);
				if (preci_flag == 0) {
					sprintf (title, "Nmovie volume_tot(A3) V(%%) <d>(A) density(m-3) N C[%s]", elem_preci[element_id]);
					fprintf (output, "%s\n", title);
//					preci_flag = 1;
				}
				V0 = VProd (region) / nMol; //dimensionless
				if (N_preci != 0) d = pow (6. * (N * V0 / N_preci) / M_PI, 1. / 3.); //dimensionless
				else d = 0.;
				densi = (double)N_preci / VProd(region); //dimensionless
				j = 0;
				DO_MOL {
					if (mol[n].preci_flag == 1 && strcmp (mol[n].elem, elem_preci[element_id]) == 0) j ++;
				}
				if (N > 0) C = (double)j / (double)N;
				else C = 0.;
				fprintf (output, "%d %e %e %e %e %d %f\n", \
					 Nmovie, N * V0 * Cube(lUnit * 1.e10), (double)N / (double)nMol * 100., d * lUnit * 1.e10, densi / Cube(lUnit), N_preci, C); //--A3, A, m-3
				printf ("%s: precipitates volume: %.3f A3, average diameter: %.3f A, density: %.3e m-3, number: %d\n", \
					elem_preci[element_id], N * V0 * Cube(lUnit * 1.e10), d * lUnit * 1.e10, densi / Cube(lUnit), N_preci); //--A3
				fclose (output);
				sprintf (filename, "../out/precipitates_%s.movie", elem_preci[element_id]);
				output = WriteFile (filename);
				if (N > 0) {
					fprintf (output, \
					 "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3:mark:I:1:C:R:1\n", \
					 N, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
				}
				DO_MOL {
					if (mol[n].preci_flag == 1) {
						fprintf (output, "%d %s %f %f %f %d %f\n", \
							 n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10), \
							 mol[n].id, mol[n].C); //--A
					}
				}
				fclose (output);
/*
				sprintf (filename, "../out/precipitates_core_%s.movie", elem_preci[element_id]);
				output = WriteFile (filename);
				if (Ncore > 0) {
					fprintf (output, \
					 "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:pos:R:3:preci_id:I:1\n", \
					 Ncore, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
				}
				for (j = 0; j < Ncore; j ++) {
					fprintf (output, "%d %f %f %f %d\n", \
						 j, coreR[j].x * (lUnit * 1.e10), coreR[j].y * (lUnit * 1.e10), coreR[j].z * (lUnit * 1.e10), \
						 core_flag[j]); //--A
				}
				fclose (output);
*/
				//identify interface
				DO_MOL {
					if (mol[n].preci_flag == 1) mol[n].interface_flag = 0;
				}
			}
/*
			//output interface precipitation file
			N = 0;
			DO_MOL {
				if (mol[n].interface_flag == 1) N ++;
			}
			sprintf (filename, "../out/precipitates_interface.dat");
			output = WriteFile (filename);
			if (preci_flag == 0) {
				fprintf (output, "Nmovie volume_tot(A3) C[%s]\n", elem_preci[0]);
				preci_flag = 1;
			}
			V0 = VProd (region) / nMol; //dimensionless
			j = 0;
			DO_MOL {
				if (mol[n].interface_flag == 1 && strcmp (mol[n].elem, elem_preci[0]) == 0) j ++;
			}
			if (N > 0) C = (double)j / (double)N;
			else C = 0;
			fprintf (output, "%d %e %f\n", Nmovie, N * V0 * Cube(lUnit * 1.e10), C); //--A3, A, m-3
			printf ("interface: precipitates volume: %.3f A3\n", \
				N * V0 * Cube(lUnit * 1.e10)); //--A3
			fclose (output);
			sprintf (filename, "../out/precipitates_interface.movie");
			output = WriteFile (filename);
			if (N > 0) {
				fprintf (output, \
				 "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3:C:R:1\n", \
				 N, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
			}
			DO_MOL {
				if (mol[n].interface_flag == 1) {
					fprintf (output, "%d %s %f %f %f %f\n", \
						 n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10), \
						 mol[n].C); //--A
				}
			}
			fclose (output);
*/
		}

		i ++;
	}

	free (coreR);
	free (core_flag);
	free (nebrTab);
	DO_MOL free (mol[n].nebr_id);
	free (mol);

	fclose (input);

	system ("rm ../out/*_SRO.dat");
	for (i = 0; i < Nelems; i ++) {
		flag_SRO = 0;
		sprintf (filename, "../out/precipitates_%s.movie", elem_preci[i]);
		CalSROMovie (filename, elemName[0], elemName[1]);
	}
	flag_SRO = 0;
	sprintf (filename, "../out/precipitates_interface.movie");
	CalSROMovie (filename, elemName[0], elemName[1]);
}

void CalSROMovie (char *file, char *element1, char *element2)
{
	FILE *input, *output;
	char filename[128], line[1024], title[1024];
	int n, N, id, temp, i, j, k, j1, j2, Nmovie, nebr_flag = 0;
	double r, tmp;
	VecR dr;;
	const int NMAX = 10000000;

	//Unit
	TUnit = 300.; //--K
	EUnit = TUnit * kB; //--J
	VUnit = EUnit / eleChar; //--eV
	lUnit = 1.e-10; //--m
	eUnit = EUnit; //--J

	input = ReadFile (file);

	//read out.movie
	i = nMol = Nmovie = stepCount = 0;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;

		// read data
		if (i == nMol) {
			Nmovie ++;
			stepCount ++;
			sscanf (line, "%d", &nMol);
			
			if (nMol > NMAX) {
				printf ("error(measurements.c): NMAX is too small\n");
				exit (1);
			}

			if (nebr_flag == 0) {
				AllocMem (Natoms, Nelems+1, int);
				AllocMem (mol, nMol, Mol);
				AllocMem (atom, nMol, Atom);
				Natoms[0] = nMol;
				DO_ATOM {
					AllocMem (atom[n].nebrLists, Ndists, NebrLists);
					for (j = 0; j < Ndists; j ++) AllocMem (atom[n].nebrLists[j].nebr, NeighbourMax, int);
				}
				AllocMem (pairDists, Ndists, double);
			} else {
				ReAllocMem (Natoms, Nelems+1, int);
				ReAllocMem (mol, nMol, Mol);
				ReAllocMem (atom, nMol, Atom);
				Natoms[0] = nMol;
				DO_ATOM {
					AllocMem (atom[n].nebrLists, Ndists, NebrLists);
					for (j = 0; j < Ndists; j ++) AllocMem (atom[n].nebrLists[j].nebr, NeighbourMax, int);
				}
				ReAllocMem (pairDists, Ndists, double);
			}

			fgets (line, 1024, input);
			sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", \
				&region.x, &tmp, &tmp, &tmp, &region.y, &tmp, &tmp, &tmp, &region.z); //--A
			if (isnormal(VProd (region)) == 0) {
				printf ("error(measurement.c): out.movie format, region(%f %f %f)\n", region.x, region.y, region.z);
				exit (1);
			}
			VScale (region, 1. / (lUnit * 1.e10)); //dimensionless
			fgets (line, 1024, input);		
			i = 0;
		}
		sscanf (line, "%d %s %lg %lg %lg", &temp, mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z); //--A
		sscanf (line, "%d %s %lg %lg %lg", &temp, atom[i].element, &atom[i].r.x, &atom[i].r.y, &atom[i].r.z);  //--A
		VScale (mol[i].r, 1. / (lUnit * 1.e10)); //dimensionless

		/*******************calcualte SRO*********************/
		if (i == nMol - 1) {
			//set initial parameter and allocate memory, build neighbour list
			rCut = 10.; // --A
			rCut /= lUnit * 1.e10; //dimensionless
			nebrTabMax = nebrTabFac * nMol;
			if (nebr_flag == 0) {
				AllocMem (nebrTab, 2 * nebrTabMax, int);
				DO_MOL AllocMem (mol[n].nebr_id, nebrTabFac, int);
			} else {
				ReAllocMem (nebrTab, 2 * nebrTabMax, int);
				DO_MOL AllocMem (mol[n].nebr_id, nebrTabFac, int);
			}
			printf ("\n-- movie: %d, build neighbour list...\n", Nmovie);
			VSCopy (cells, 1. / (rCut + rNebrShell), region);
			AllocMem (cellList, VProd (cells) + nMol, int);
			BuildNebrList ();
			free (cellList);
			//atom numbers of different elements
			printf ("atom number: ");
			for (j = 1; j <= Nelems; j ++) {
				Natoms[j] = 0;
				DO_ATOM {
					if (strcmp(atom[n].element, elemName[j-1]) == 0) Natoms[j] ++;
				}
				printf ("%s %d %.2f%% ", elemName[j-1], Natoms[j], (double)(Natoms[j])/(double)(Natoms[0])*100.);
			}
			printf ("\n");
			Get_pair_distance ();
			//get i-th sphere neighbour atoms and pairs lists
			Get_neighbour ();
			nebr_flag = 1;
			
			//calcualte short rang order (SRO)
			sprintf (filename, "../out/%s_SRO.dat", file);
			ShortRangeOrder (element1, element2, filename);
			printf ("\n");

			//free the points in the structure
			DO_MOL free (mol[n].nebr_id);
			DO_ATOM {
				for (j = 0; j < Ndists; j ++) free (atom[n].nebrLists[j].nebr);
				free (atom[n].nebrLists);
			}
		}
		i ++;
	}

	fclose (input);
}

void CalSROMovie_TOTAL ()
{
	FILE *input, *output;
	char line[1024], filename[128];
	int N, i, j, k, flag = 0, id, mem_flag = 0;
	double dr_threshold, dr, temp, drmin1, drmin2;
	Mol *mol;
	VecR region, Vdr;

	int *cellList;
	VecI cells;
	VecR invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr, rNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

	int m, nebr, alpha[20];
	int NiAB, NiA, NA;
	double alphai, piAB, cB, PiAB, ci;

	double ZFe, Z, alpha1, ZFe1, Z1, alpha2, ZFe2, Z2, cCr, beta;

	sprintf (filename, "../out/md/out.movie");
	input = fopen (filename, "r");

//	printf ("the radius of shell for SRO is:(eg. 3.)\n");
//	scanf ("%lg", &dr_threshold);
	dr_threshold = 6.;

	rNebr = dr_threshold;
	stepCount = 0;
	while (1) {
		if (fgets (line, 1024, input) == 0) break;
		sscanf (line, "%d", &N);
		nMol = N;
		fgets (line, 1024, input);
		sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", &region.x, &temp, &temp, &temp, &region.y, &temp, &temp, &temp, &region.z);
		if (mem_flag == 0) {
			printf ("region:[%f %f %f]\n", region.x, region.y, region.z); 
			printf ("building neighbour list...\n");
			AllocMem (mol, N, Mol);
			for (i = 0; i < N; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z);
			}
			k = 0;
			for (i = 0; i < N; i ++) mol[i].Nnebr = 0;
			nebrTabMax = 200 * nMol;
			VSCopy (cells, 1. / rNebr, region);
			printf ("cells: [%d %d %d]\n", cells.x, cells.y, cells.z);
			AllocMem (cellList, VProd (cells) + nMol, int);
			AllocMem (nebrTab, 2 * nebrTabMax, int);
			rrNebr = Sqr (rNebr);
			VDiv (invWid, cells, region);
			for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
			DO_MOL{
		//		VSAdd (rs, mol[n].r, 0.5, region);
				VCopy (rs, mol[n].r);
				VMul (cc, rs, invWid);
				c = VLinear (cc, cells) + nMol;
				if (c >= (VProd(cells) + nMol) || c < 0) {
					printf ("error(nebrList.c): cellList is not big enough(measurement.c),\n \
						cc:(%d %d %d), cells(%d %d %d)\n \
						rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
						rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
						region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
					PrintMovie ();
					exit (1);
				}
				cellList[n] = cellList[c];
				cellList[c] = n;
			}
			nebrTabLen = 0;
			for (m1z = 0; m1z < cells.z; m1z ++){
				for (m1y = 0; m1y < cells.y; m1y ++){
					for (m1x = 0; m1x < cells.x; m1x ++){
						VSet (m1v, m1x, m1y, m1z);
						m1 = VLinear (m1v, cells) + nMol;
						for (offset = 0; offset < N_OFFSET; offset ++){
							VAdd (m2v, m1v, vOff[offset]);
							VZero (shift);
							VCellWrapAll ();
							m2 = VLinear (m2v, cells) + nMol;
							DO_CELL (j1, m1){
								DO_CELL (j2,m2){
									if (m1 != m2 || j2 < j1){
										VSub (Vdr, mol[j1].r, mol[j2].r);
										VVSub (Vdr, shift);
										if (VLenSq (Vdr) < rrNebr){
											if (nebrTabLen >= nebrTabMax) {
												printf ("error(measurement.c): ERR_TOO_MANY_NEBRS");
												exit (1);
											}
											nebrTab[2 * nebrTabLen] = j1;
											nebrTab[2 * nebrTabLen + 1] = j2;
											++ nebrTabLen;
											i = mol[j1].Nnebr;
											j = mol[j2].Nnebr;
											mol[j1].nebr[i] = j2;
											mol[j2].nebr[j] = j1;
											mol[j1].Nnebr ++;
											mol[j2].Nnebr ++;
										}
									}
								}
							}
						}
					}
				}
			}
			//find rmin1 and rmin2
			drmin1 = drmin2 = 1.e10;
			for (i = 0; i < nebrTabLen; i ++) {
				j1 = nebrTab[2 * i];
				j2 = nebrTab[2 * i + 1];
				VSub (Vdr, mol[j1].r, mol[j2].r);
				PeriodBoundary (Vdr);
				dr = VLen (Vdr);
				if (drmin1 > dr) drmin1 = dr;
			}
			for (i = 0; i < nebrTabLen; i ++) {
				j1 = nebrTab[2 * i];
				j2 = nebrTab[2 * i + 1];
				VSub (Vdr, mol[j1].r, mol[j2].r);
				PeriodBoundary (Vdr);
				dr = VLen (Vdr);
				if ((drmin2 > dr) && (dr - drmin1 > 0.01)) drmin2 = dr;
			}
			printf ("the nearest proximity distance is %f\nthe secondary proximity distance is %f\n", drmin1, drmin2);

			mem_flag = 1;
		} else {
			for (i = 0; i < N; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &temp, &temp, &temp);
			}
		}

		//average SRO
		// alpha = 1 - ZFe / (Z * (1 - cCr))
		// beta = (8 * alpha1 + 6 * alpha2) / 14
		ZFe = Z = 0;
		for (n = 0; n < N; n ++) {
			if (strcmp(mol[n].elem, elemName[0]) != 0) ZFe ++;
			Z ++;
		}
		cCr = 1. - (double)ZFe / (double)Z;
		beta = 0.;
		for (m = -10; m < 10; m ++) alpha[m] = 0;
		for (n = 0; n < N; n ++) {
			if (strcmp(mol[n].elem, elemName[0]) == 0) {
				//the  nearest shell
				ZFe1 = Z1 = ZFe2 = Z2 = 0;
				for (m = 0; m < mol[n].Nnebr; m++) {
					nebr = mol[n].nebr[m];
					VSub (Vdr, mol[n].r, mol[nebr].r);
					PeriodBoundary (Vdr);
					dr = VLen (Vdr);
					if (fabs(dr - drmin1) < 0.01) {
						Z1 ++;
						if (strcmp(mol[nebr].elem, elemName[0]) != 0) ZFe1 ++;
					}
					if (fabs(dr - drmin2) < 0.01) {
						Z2 ++;
						if (strcmp(mol[nebr].elem, elemName[0]) != 0) ZFe2 ++;
					}
				}
				alpha1 = 1. - ZFe1 / (Z1 * (1. - cCr));
				alpha2 = 1. - ZFe2 / (Z2 * (1. - cCr));
				alphai =  (8. * alpha1 + 6. * alpha2) / 14.;
				for (m = -10; m < 10; m ++) {
					if (alphai > (double)m * 0.1 && alphai < (double)(m + 1) * 0.1) alpha[m] ++;
				}
				beta += alphai;
			} 
		}
		beta /= (double)(Z - ZFe);
		printf ("%lld SRO: %f\n", stepCount, beta);

		sprintf (filename, "../out/sro.dat");
		output = WriteFile (filename);
		if (flag_SRO == 0) fprintf (output, "step SRO\n");
		fprintf (output, "%lld %f\n", stepCount, beta);
		fclose (output);
		output = fopen ("../out/sro_profile.dat", "a+");
		fprintf (output, "SRO N Fraction_of_atoms\n");
		for (m = -10; m < 10; m ++) fprintf (output, "%f %d %f\n", ((double)m + 0.5) * 0.1, alpha[m], (double)alpha[m] / (double)(Z - ZFe) * 100.);
		fclose (output);
/*
		// alphai = 1 - PiAB / cB, cB is the concentration of atoms of B sort.
		// PiAB = NiAB / NiA, PiAB is the probability to find a B atom on the i-th sphere of coordination of an A atom located at the origin of coordinates.
		NiAB = 0;
		NiA = 0;
		NA = 0;
		i = 0;
		for (n = 0; n < N; n ++) {
			if (strcmp(mol[n].elem, "Cr") == 0) {
				NA ++;
				NiA += mol[n].Nnebr;
				for (m = 0; m < mol[n].Nnebr; m++) {
					nebr = mol[n].nebr[m];
					if (strcmp(mol[nebr].elem, "Fe") == 0) NiAB ++;
				}
			} 
		}
		PiAB = (double)NiAB / (double)NiA;
		cB = 1. - (double)NA / (double)N;
		alphai = 1. - PiAB / cB;

		printf ("%lld SRO: %f\n", stepCount, alphai);

		sprintf (filename, "../out/sro.dat");
		output = WriteFile (filename);
		if (flag_SRO == 0) fprintf (output, "step SRO\n");
		fprintf (output, "%lld %f\n", stepCount, alphai);
		fclose (output);

		//SRO profile
		for (m = -10; m < 10; m ++) alpha[m] = 0;
		NA = 0;
		for (n = 0; n < N; n ++) {
			if (strcmp(mol[n].elem, "Cr") == 0) {
				NA ++;
				NiAB = 0;
				for (m = 0; m < mol[n].Nnebr; m++) {
					nebr = mol[n].nebr[m];
					if (strcmp(mol[nebr].elem, "Fe") == 0) NiAB ++;
				}
				PiAB = (double)NiAB / (double)(mol[n].Nnebr);
				alphai = 1. - PiAB / cB;
				for (m = -10; m < 10; m ++) {
					if (alphai > (double)m * 0.1 && alphai < (double)(m + 1) * 0.1) alpha[m] ++;
				}
//				printf ("%f\n", alphai);
			} 
		}
		output = fopen ("../out/sro_profile.dat", "a+");
		fprintf (output, "SRO N Fraction_of_atoms\n");
		for (m = -10; m < 10; m ++) fprintf (output, "%f %d %f\n", ((double)m + 0.5) * 0.1, alpha[m], (double)alpha[m] / (double)NA * 100.);
		fclose (output);
*/

		flag_SRO = 1;
		stepCount ++;
	}

	free (cellList);
	free (nebrTab);
	free (mol);
	fclose (input);
}

void CalAutocorrelationFunction ()
{
	FILE *input, *output;
	char line[1024], filename[128];
	int i, j, k, flag = 0, id, mem_flag = 0;
	double dr, temp;
	Mol *mol;
	VecR region, Vdr;

	int *cellList;
	VecI cells;
	VecR invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr, rNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

	int m, NA, NAK, N, NK, Ncenter;
	double rmax, rkmax, r, rk, Cr, Crk, C0, Rk, *Rkavg, sigma2;
	VecR center;

	sprintf (filename, "../out/md/out.movie");
	input = fopen (filename, "r");
	sprintf (filename, "../out/acf.dat");

	printf ("the radius of max shell for ACF is (eg. 3.) [A]:\n");
	scanf ("%lg", &rmax);
	printf ("the number of center points is:\n");
	scanf ("%d", &Ncenter);

	stepCount = 0;
	while (1) {
		if (fgets (line, 1024, input) == 0) break;
		sscanf (line, "%d", &nMol);
		fgets (line, 1024, input);
		sscanf (line, "Lattice=\"%lg 0.0 0.0 0.0 %lg 0.0 0.0 0.0 %lg", &region.x, &region.y, &region.z);
		if (mem_flag == 0) {
			AllocMem (mol, nMol, Mol);
			AllocMem (Rkavg, (int)rmax+1, double);
			for (i = 0; i < nMol; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z);
			}
			mem_flag = 1;
		} else {
			for (i = 0; i < nMol; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &temp, &temp, &temp);
			}
		}

		//calculate ACF
		NA = 0;
		for (i = 0; i < nMol; i ++) {
			if (strcmp (mol[i].elem, elemName[0]) == 0) NA ++;
		}
		C0 = (double)NA / (double)nMol;
		printf ("C0 = %f\n", C0);
		srand ((unsigned int)time(0));
		for (k = 0; k <= (int)rmax; k ++) Rkavg[k] = 0.;
		for (m = 0; m < Ncenter; m ++) {
			center.x = (double)rand () / (double)RAND_MAX * region.x;
			center.y = (double)rand () / (double)RAND_MAX * region.y;
			center.z = (double)rand () / (double)RAND_MAX * region.z;
			printf ("%d: center: [%f %f %f]\n", m, center.x, center.y, center.z);
			output = fopen (filename, "a+");
			for (k = 0; k <= (int)rmax; k ++) {
				rkmax = rmax - k;
				Rk = 0.;
				for (j = 0; j <= (int)rkmax; j ++) {
					r = j;
					rk = r + k;
					NA = NAK = N = NK = 0;
					#pragma omp parallel for reduction (+:N, NA, NK, NAK) private (id, Vdr, dr) num_threads (nthreads)
					for (i = 0; i < nMol; i ++) {
						id = i;
						VSub (Vdr, center, mol[id].r);
						PeriodBoundary (Vdr);
						dr = VLen (Vdr);
						if (dr < r && dr > r - 1) {
							N ++;
							if (strcmp (mol[id].elem, elemName[0]) == 0) NA ++;
						}
						if (dr < rk && dr > rk - 1) {
							NK ++;
							if (strcmp (mol[id].elem, elemName[0]) == 0) NAK ++; 
						}
					}
					if (N == 0 || NK == 0) continue;
					Cr = (double)NA / (double)N;
					Crk = (double)NAK / (double)NK;
					Rk += (Cr - C0) * (Crk - C0);
				}
				sigma2 = 0.;
				for (j = 0; j <= (int)rmax; j ++) {
					r = j;
					NA = N = 0;
					#pragma omp parallel for reduction (+:N, NA) private (id, Vdr, dr) num_threads (nthreads)
					for (i = 0; i < nMol; i ++) {
						id = i;
						VSub (Vdr, center, mol[id].r);
						PeriodBoundary (Vdr);
						dr = VLen (Vdr);
						if (dr < r && dr > r - 1) {
							N ++;
							if (strcmp (mol[id].elem, elemName[0]) == 0) NA ++;
						}
					}
					if (N == 0) continue;
					Cr = (double)NA / (double)N;
					sigma2 += Sqr (Cr - C0);
				}
				Rk = Rk / sigma2;
				Rkavg[k] += Rk / (double)Ncenter;
//				printf ("%d %f\n", k, Rk);
			}
		}
		fprintf (output, "Distance(nm) Rk\n");
		for (k = 0; k <= (int)rmax; k ++) {
			printf ("%f %f\n", k / 10., Rkavg[k]);
			fprintf (output, "%f %f\n", k / 10., Rkavg[k]);
		}
		fclose (output);

		printf ("\nstep: %lld\n", stepCount);
		stepCount ++;
	}

	free (mol);
	free (Rkavg);
	fclose (input);
}

void ClusterExtract ()
{
	FILE *input, *output, *output1;
	char line[1024], filename[128];
	int N, i, j, k, l, flag = 0, id, mem_flag = 0, Nca, *ca;
	double dr_threshold, dr, temp, drmin1, drmin2, drmin3;
	VecR Vdr;

	int *cellList;
	VecI cells;
	VecR invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr, rNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, m, nebr, offset;

	int Nc, sym, cluster_type, NCr, NFe;
	Mol *ce;
	VecR move, move1, cenb;

	sprintf (filename, "../out/md/out.movie");
	input = fopen (filename, "r");

	stepCount = 0;
	while (1) {
		if (fgets (line, 1024, input) == 0) break;
		sscanf (line, "%d", &N);
		nMol = N;
		fgets (line, 1024, input);
		sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", &region.x, &temp, &temp, &temp, &region.y, &temp, &temp, &temp, &region.z);
		if (mem_flag == 0) {
			printf ("region:[%f %f %f]\n", region.x, region.y, region.z); 
			printf ("building neighbour list...\n");
			AllocMem (mol, N, Mol);
			for (i = 0; i < N; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &mol[i].r.x, &mol[i].r.y, &mol[i].r.z);
				if (strcmp(mol[i].elem, "vacancy") == 0) sprintf (mol[i].elem, "%s", elemName[0]);
			}

			//find rmin1, rmin2 and rmin3
			drmin1 = drmin2 = drmin3 = 1.e10;
			for (i = 1; i < nMol; i ++) {
				j1 = 0;
				j2 = i;
				VSub (Vdr, mol[j1].r, mol[j2].r);
				PeriodBoundary (Vdr);
				dr = VLen (Vdr);
				if (drmin1 > dr) drmin1 = dr;
			}
			for (i = 1; i < nMol; i ++) {
				j1 = 0;
				j2 = i;
				VSub (Vdr, mol[j1].r, mol[j2].r);
				PeriodBoundary (Vdr);
				dr = VLen (Vdr);
				if ((drmin2 > dr) && (dr - drmin1 > 0.01)) drmin2 = dr;
			}
			for (i = 1; i < nMol; i ++) {
				j1 = 0;
				j2 = i;
				VSub (Vdr, mol[j1].r, mol[j2].r);
				PeriodBoundary (Vdr);
				dr = VLen (Vdr);
				if ((drmin3 > dr) && (dr - drmin2 > 0.01)) drmin3 = dr;
			}
			printf ("the nearest proximity distance is %f\nthe secondary proximity distance is %f\n", drmin1, drmin2);
			printf ("the third proximity distance is %f\n", drmin3);

			//build neighbour list
			printf ("\nchoose cluster type: \n[\n 1(the nearest shell);\n 2(the second nearest shell);\n 3(the third nearest shell)\n]:\n");
//			scanf ("%d", &cluster_type);
			cluster_type = (int)ith_nebrR[0];
			printf ("cluster_type: %d\n", cluster_type);
			switch (cluster_type) {
				case 0: rNebr = 3.464 * drmin2; break;
				case 1: rNebr = drmin1; break;
				case 2: rNebr = drmin2; break;
				case 3: rNebr = drmin3; break;
				default: break;
			}
			rNebr += 0.01;
			for (i = 0; i < N; i ++) mol[i].Nnebr = 0;
			nebrTabMax = 200 * nMol;
			VSCopy (cells, 1. / rNebr, region);
			printf ("cells: [%d %d %d]\n", cells.x, cells.y, cells.z);
			AllocMem (cellList, VProd (cells) + nMol, int);
			AllocMem (nebrTab, 2 * nebrTabMax, int);
			rrNebr = Sqr (rNebr);
			VDiv (invWid, cells, region);
			for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
			DO_MOL{
//				VSAdd (rs, mol[n].r, 0.5, region);
				VCopy (rs, mol[n].r);
				VMul (cc, rs, invWid);
				c = VLinear (cc, cells) + nMol;
				if (c >= (VProd(cells) + nMol) || c < 0) {
					printf ("error(nebrList.c): cellList is not big enough(measurement.c),\n \
						cc:(%d %d %d), cells(%d %d %d)\n \
						rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
						rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
						region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
					PrintMovie ();
					exit (1);
				}
				cellList[n] = cellList[c];
				cellList[c] = n;
			}
			nebrTabLen = 0;
			for (m1z = 0; m1z < cells.z; m1z ++){
				for (m1y = 0; m1y < cells.y; m1y ++){
					for (m1x = 0; m1x < cells.x; m1x ++){
						VSet (m1v, m1x, m1y, m1z);
						m1 = VLinear (m1v, cells) + nMol;
						for (offset = 0; offset < N_OFFSET; offset ++){
							VAdd (m2v, m1v, vOff[offset]);
							VZero (shift);
							VCellWrapAll ();
							m2 = VLinear (m2v, cells) + nMol;
							DO_CELL (j1, m1){
								DO_CELL (j2,m2){
									if (m1 != m2 || j2 < j1){
										VSub (Vdr, mol[j1].r, mol[j2].r);
										VVSub (Vdr, shift);
										if (VLenSq (Vdr) < rrNebr){
											if (nebrTabLen >= nebrTabMax) {
												printf ("error(measurement.c): ERR_TOO_MANY_NEBRS");
												exit (1);
											}
											nebrTab[2 * nebrTabLen] = j1;
											nebrTab[2 * nebrTabLen + 1] = j2;
											++ nebrTabLen;
											i = mol[j1].Nnebr;
											j = mol[j2].Nnebr;
											if (cluster_type != 0) {
												mol[j1].nebr[i] = j2;
												mol[j2].nebr[j] = j1;
												mol[j1].Nnebr ++;
												mol[j2].Nnebr ++;
											} else {
												if ((Vdr.x > -0.01 && Vdr.x < 2. * drmin2 - 0.01) &&
												    (Vdr.y > -0.01 && Vdr.y < 2. * drmin2 - 0.01) &&
												    (Vdr.z > -0.01 && Vdr.z < 2. * drmin2 - 0.01)) {
													mol[j2].nebr[j] = j1;
													mol[j2].Nnebr ++;
												} else if ((Vdr.x < 0.01 && Vdr.x > -2 * drmin2 + 0.01) &&
												    (Vdr.y < 0.01 && Vdr.y > -2. * drmin2 + 0.01) &&
												    (Vdr.z < 0.01 && Vdr.z > -2. * drmin2 + 0.01)) {
													mol[j1].nebr[i] = j2;
													mol[j1].Nnebr ++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			mem_flag = 1;
		} else {
			for (i = 0; i < N; i ++) {
				fgets (line, 1024, input);
				sscanf (line, "%lg %s %lg %lg %lg", &temp, mol[i].elem, &temp, &temp, &temp);
			}
		}

		Nc = 0;
		AllocMem (ce, 1, Mol);
		Nca = mol[0].Nnebr + 2;
		for (n = 0; n < N; n ++) {
			if (mol[n].Nnebr + 2 != Nca) {
				printf ("error: box region may be wrong\n");
				exit (1);
			}
			i = 0;
			if (strcmp(mol[n].elem, elem_preci[0]) == 0) i ++;
			for (m = 0; m < mol[n].Nnebr; m++) {
				nebr = mol[n].nebr[m];
				if (strcmp(mol[nebr].elem, elem_preci[0]) == 0) i ++;
			}
			mol[n].C = (double)i / (double)(mol[n].Nnebr + 1);
			flag = 1; //if adding cluster to lib
			for (i = 0; i < Nc; i ++) {
				if (fabs(mol[n].C - ce[i].C) < 1.e-100) {
					sym = Symmetry (mol[n], ce[i]);
					if (sym == 1) {
						flag = 0;
						ce[i].Nsym ++;
					}
				}
			}
			if (flag == 1) {
				Nc ++;
				ReAllocMem (ce, Nc, Mol);
				ce[Nc-1] = mol[n];
				ce[Nc-1].Nsym = 1;
//				printf ("%d %.3f\n", Nc-1, ce[Nc-1].C);
			}
			if (n % (N / 100) == 0) printf ("clusters: %d, %.0f%%\n", Nc, n * 100. / N);
		}

		AllocMem (ca, Nca, int);
		for (i = 0; i < Nca; i ++) {
			ca[i] = 0;
		}
		for (n = 0; n < N; n ++) {
			id = mol[n].C * (double)(Nca - 1);
			ca[id] ++;
		}
		output = fopen ("../out/CE/ce.dat", "w");
		printf ("\nC Nce proportion\n");
		fprintf (output, "C Nce proportion\n");
		for (i = 0; i < Nca; i ++) {
			printf ("%.3f %d %.3f\n", (double)i / (double)(Nca - 1), ca[i], (double)ca[i] / (double)N);
			fprintf (output, "%.3f %d %.3f\n", (double)i / (double)(Nca - 1), ca[i], (double)ca[i] / (double)N);
		}
		free (ca);
		fclose (output);

		for (i = 0; i < Nc; i ++) {
			sprintf (filename, "../out/CE/ce_%.3f.movie", ce[i].C);
			output = fopen (filename, "a+");
			fprintf (output, "%d\nProperties=id:I:1:species:S:1:pos:R:3\n", ce[i].Nnebr+1);
			sprintf (filename, "../out/CE/cell_%.3f.movie", ce[i].C);
			output1 = fopen (filename, "a+");
			fprintf (output1, "%d\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" Properties=id:I:1:species:S:1:pos:R:3\n", ce[i].Nnebr+2, drmin2 * 2., drmin2 * 2., drmin2 * 2.);
			VSCopy (move, -1., ce[i].r);
			fprintf (output, "%d %s %e %e %e\n", 0, ce[i].elem, 0., 0., 0.);
			fprintf (output1, "%d %s %e %e %e\n", 0, ce[i].elem, drmin2, drmin2, drmin2);
			j = 1;
			NCr = NFe = 0;
			
			for (m = 0; m < ce[i].Nnebr; m ++) {
				//cluster exterct
				nebr = ce[i].nebr[m];
				VAdd (cenb, mol[nebr].r, move);
				if (cenb.x > region.x * 0.5) cenb.x -= region.x;
				else if (cenb.x < -region.x * 0.5) cenb.x += region.x;
				if (cenb.y > region.y * 0.5) cenb.y -= region.y;
				else if (cenb.y < -region.y * 0.5) cenb.y += region.y;
				if (cenb.z > region.z * 0.5) cenb.z -= region.z;
				else if (cenb.z < -region.z * 0.5) cenb.z += region.z;
				fprintf (output, "%d %s %e %e %e\n", m+1, mol[nebr].elem, cenb.x, cenb.y, cenb.z);

				//cell exterct
				VSet (move1, drmin2, drmin2, drmin2);
				VVAdd (cenb, move1);
				if (fabs(cenb.x - 2. * drmin2) < 0.01 || fabs(cenb.y - 2. * drmin2) < 0.01  || fabs(cenb.z - 2. * drmin2) < 0.01) {
					if (strcmp (mol[nebr].elem, elem_preci[0]) == 0) NCr ++;
					else NFe ++;
					continue;
				}
				fprintf (output1, "%d %s %e %e %e\n", j, mol[nebr].elem, cenb.x, cenb.y, cenb.z);
				j ++;
			}
			
			if (ce[i].C > 0.5) fprintf (output1, "%d %s %e %e %e\n", j, elem_preci[0], 0., 0., 0.);
			else fprintf (output1, "%d Fe %e %e %e\n", j, 0., 0., 0.);
			if (NCr >= 1) fprintf (output1, "%d Cr %e %e %e\n", j + 1, 0., 0., drmin2);
			else fprintf (output1, "%d Fe %e %e %e\n", j + 1, 0., 0., drmin2);
			if (NCr >= 2) fprintf (output1, "%d Cr %e %e %e\n", j + 2, 0., drmin2, 0.);
			else fprintf (output1, "%d Fe %e %e %e\n", j + 2, 0., drmin2, 0.);
			if (NCr >= 3) fprintf (output1, "%d Cr %e %e %e\n", j + 3, drmin2, 0., 0.);
			else fprintf (output1, "%d Fe %e %e %e\n", j + 3, drmin2, 0., 0.);

			fclose (output);
			fclose (output1);
			sprintf (filename, "../out/CE/ce_%.3f.dat", ce[i].C);
			output = fopen (filename, "a+");
			fprintf (output, "%d\n", ce[i].Nsym);
			fclose (output);
		}
		free (ce);


		flag_SRO = 1;
		stepCount ++;
	}

	free (cellList);
	free (nebrTab);
	free (mol);
	fclose (input);
}

int Symmetry (Mol mol1, Mol mol2)
{
	int flag, i, j, nebr1, nebr2, sym;
	VecR Vdr1, Vdr2, Vdr;
	double dr, temp;

	if (strcmp (mol1.elem, mol2.elem) != 0) return (0);

	sym = 0;

	//[1 0 0]
	//[0 1 0]
	//[0 0 1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[1 0 0]
	//[0 1 0]
	//[0 0 -1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.z *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[1 0 0]
	//[0 -1 0]
	//[0 0 1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.y *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[-1 0 0]
	//[0 1 0]
	//[0 0 1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.x *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[-1 0 0]
	//[0 -1 0]
	//[0 0 -1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.x *= -1.;
		Vdr1.y *= -1.;
		Vdr1.z *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[1 0 0]
	//[0 -1 0]
	//[0 0 -1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.y *= -1.;
		Vdr1.z *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[-1 0 0]
	//[0 1 0]
	//[0 0 -1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.x *= -1.;
		Vdr1.z *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[-1 0 0]
	//[0 -1 0]
	//[0 0 1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		Vdr1.x *= -1.;
		Vdr1.y *= -1.;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;
/*
	//[1 0 0]
	//[0 0 1]
	//[0 1 0]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		temp = Vdr1.y;
		Vdr1.y = Vdr1.z;
		Vdr1.z = temp;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[0 0 1]
	//[0 1 0]
	//[1 0 0]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		temp = Vdr1.x;
		Vdr1.x = Vdr1.z;
		Vdr1.z = temp;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;

	//[0 1 0]
	//[1 0 0]
	//[0 0 1]
	flag = 1;
	for (i = 0; i < mol1.Nnebr; i ++) {
		nebr1 = mol1.nebr[i];
		VSub (Vdr1, mol1.r, mol[nebr1].r);
		PeriodBoundary (Vdr1);
		temp = Vdr1.x;
		Vdr1.x = Vdr1.y;
		Vdr1.y = temp;
		for (j = 0; j < mol2.Nnebr; j ++) {
			nebr2 = mol2.nebr[j];
			VSub (Vdr2, mol2.r, mol[nebr2].r);
			PeriodBoundary (Vdr2);
			VSub (Vdr, Vdr1, Vdr2);
			dr = VLen (Vdr);
			if (fabs(dr) < 0.01)
				if (strcmp (mol[nebr1].elem, mol[nebr2].elem) != 0) flag = 0;
		}
	}
	if (flag == 1) sym = 1;
*/
	return (sym);
}

void KMC_Cluster ()
{
	int i, j, k, NCr, N, id, Nca, *ca, Nca_goal, loop, loop_max, Nelem, Ngoal, count, Ndr, *Nelem_atom, Norder;
	double drmin1, drmin2, drmin3, dr, C, **cap, **cap_goal, temp, diff, diff_old, prob, *Celem;
	VecR Vdr, region_goal;
	FILE *input, *output;
	char line[1024], elem[8], **element, filename[128];

	int *cellList;
	VecI cells;
	VecR invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr, rNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, m, nebr, offset, error;

	Mol *goal;

	//obtain goal strcture parameters
/*
	printf ("input elements number:\n");
	scanf ("%d", &Nelem);
	AllocMem2 (element, Nelem, 8, char);
	AllocMem (Celem, Nelem, double);
	printf ("input elements name:\n");
	for (i = 0; i < Nelem; i ++) {
		scanf ("%s", element[i]);
	}
	printf ("input the loop number:(eg. 1000000)\n");
	scanf ("%d", &loop_max);
	printf ("input the probilibity of exchanging atoms if diff > 0:(eg. 0.00001)\n");
	scanf ("%lg", &prob);
	printf ("input the error type: standard deviation(1) or absolute difference(2)\n");
	scanf ("%d", &error);
	printf ("input the number of shell for the cluster radius:(1, 2 or 3):\n");
	scanf ("%d", &Ndr);
	printf ("inut whether the goal atom structure is order(1) or disorder (0):\n");
	scanf ("%d", &Norder);
*/
	Nelem = Nelems;
	AllocMem2 (element, Nelem, 8, char);
	AllocMem (Celem, Nelem, double);
	for (i = 0; i < Nelem; i ++) {
		sprintf (element[i], "%s", elemName[i]);
	}
	loop_max = NLOOP;
	Ndr = (int)ith_nebrR[0];
	if (strcmp(cell_order, "disorder") == 0) Norder = 0;
	else Norder = 1;
	prob = 0.00001;
	error = error_map;

	printf ("\ngoal structure:\n");
	sprintf (filename, "../in/%s", filename_map);
	input = fopen (filename, "r");
	fgets (line, 1024, input);
	sscanf (line, "%d", &Ngoal);
	AllocMem (goal, Ngoal, Mol);
	for (n = 0; n < Ngoal; n ++) AllocMem (goal[n].Celem, Nelem, double);
	fgets (line, 1024, input);
	sscanf (line, "Lattice=\"%lg 0.0 0.0 0.0 %lg 0.0 0.0 0.0 %lg", &region_goal.x, &region_goal.y, &region_goal.z);
	printf ("goal atoms: %d, box: [%f %f %f]\n", Ngoal, region_goal.x, region_goal.y, region_goal.z);
	for (n = 0; n < Ngoal; n ++) {
		fgets (line, 1024, input);
		sscanf (line, "%s %lg %lg %lg", goal[n].elem, &goal[n].r.x, &goal[n].r.y, &goal[n].r.z);
	}
	fclose (input);
	for (i = 0; i < Nelem; i ++) {
		count = 0;
		for (n = 0; n < Ngoal; n ++) {
			if (strcmp (goal[n].elem, element[i]) == 0) count ++;
		}
		Celem[i] = (double)count / (double)Ngoal;
		printf ("%s: %.3f\n", element[i], Celem[i]);
	}
	//find rmin1, rmin2 and rmin3
	drmin1 = drmin2 = drmin3 = 1.e10;
	for (i = 1; i < Ngoal; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, goal[j1].r, goal[j2].r);
		PeriodBoundaryGoal (Vdr);
		dr = VLen (Vdr);
		if (drmin1 > dr) drmin1 = dr;
	}
	for (i = 1; i < Ngoal; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, goal[j1].r, goal[j2].r);
		PeriodBoundaryGoal (Vdr);
		dr = VLen (Vdr);
		if ((drmin2 > dr) && (dr - drmin1 > 0.01)) drmin2 = dr;
	}
	for (i = 1; i < Ngoal; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, goal[j1].r, goal[j2].r);
		PeriodBoundaryGoal (Vdr);
		dr = VLen (Vdr);
		if ((drmin3 > dr) && (dr - drmin2 > 0.01)) drmin3 = dr;
	}
	printf ("the nearest proximity distance is %f\nthe secondary proximity distance is %f\n", drmin1, drmin2);
	printf ("the third proximity distance is %f\n", drmin3);
			//build cell list
			if (Ndr == 1) rNebr = drmin1 + 0.01;
			else if (Ndr == 2) rNebr = drmin2 + 0.01;
			else rNebr = drmin3 + 0.01;
			for (i = 0; i < Ngoal; i ++) goal[i].Nnebr = 0;
			nebrTabMax = 200 * Ngoal;
			VSCopy (cells, 1. / rNebr, region_goal);
			printf ("cells: [%d %d %d]\n", cells.x, cells.y, cells.z);
			AllocMem (cellList, VProd (cells) + Ngoal, int);
			AllocMem (nebrTab, 2 * nebrTabMax, int);
			rrNebr = Sqr (rNebr);
			VDiv (invWid, cells, region_goal);
			for (n = Ngoal; n < Ngoal + VProd (cells); n ++) cellList[n] = -1;
			for (n = 0; n < Ngoal; n ++) {
//				VSAdd (rs, mol[n].r, 0.5, region);
				VCopy (rs, goal[n].r);
				VMul (cc, rs, invWid);
				c = VLinear (cc, cells) + Ngoal;
				if (c >= (VProd(cells) + Ngoal) || c < 0) {
					printf ("error(nebrList.c): cellList is not big enough(measurement.c),\n \
						cc:(%d %d %d), cells(%d %d %d)\n \
						rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
						rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
						region_goal.x * lUnit * 1.e10, region_goal.y * lUnit * 1.e10, region_goal.z * lUnit * 1.e10);
					exit (1);
				}
				cellList[n] = cellList[c];
				cellList[c] = n;
			}
			nebrTabLen = 0;
			for (m1z = 0; m1z < cells.z; m1z ++){
				for (m1y = 0; m1y < cells.y; m1y ++){
					for (m1x = 0; m1x < cells.x; m1x ++){
						VSet (m1v, m1x, m1y, m1z);
						m1 = VLinear (m1v, cells) + Ngoal;
						for (offset = 0; offset < N_OFFSET; offset ++){
							VAdd (m2v, m1v, vOff[offset]);
							VZero (shift);
							VCellWrapAllGoal ();
							m2 = VLinear (m2v, cells) + Ngoal;
							DO_CELL (j1, m1){
								DO_CELL (j2,m2){
									if (m1 != m2 || j2 < j1){
										VSub (Vdr, goal[j1].r, goal[j2].r);
										VVSub (Vdr, shift);
										if (VLenSq (Vdr) < rrNebr){
											if (nebrTabLen >= nebrTabMax) {
												printf ("error(measurement.c): ERR_TOO_MANY_NEBRS");
												exit (1);
											}
											nebrTab[2 * nebrTabLen] = j1;
											nebrTab[2 * nebrTabLen + 1] = j2;
											++ nebrTabLen;
											i = goal[j1].Nnebr;
											j = goal[j2].Nnebr;
											goal[j1].nebr[i] = j2;
											goal[j2].nebr[j] = j1;
											goal[j1].Nnebr ++;
											goal[j2].Nnebr ++;
										}
									}
								}
							}
						}
					}
				}
			}
			free (cellList);
			free (nebrTab);
	for (j = 0; j < Nelem; j ++) {
		for (n = 0; n < Ngoal; n ++) {
			goal[n].Celem[j] = 0.;
			if (strcmp(goal[n].elem, element[j]) == 0) goal[n].Celem[j] += 1.;
			for (i = 0; i < goal[n].Nnebr; i ++) {
				id = goal[n].nebr[i];
				if (strcmp(goal[id].elem, element[j]) == 0)  goal[n].Celem[j] += 1.;
			}
			goal[n].Celem[j] /= (goal[n].Nnebr + 1.);
		}
	}
	/*
	for (n = 0; n < Ngoal; n ++) {
		for (j = 0; j < Nelem; j ++) {
			printf ("%s %.3f ", element[j], goal[n].Celem[j]);
		}
		printf ("\n");
	}
	*/
	Nca_goal = goal[0].Nnebr + 2;
	AllocMem2 (cap_goal, Nelem, Nca_goal, double); //the proportion of each kind of clusters
	
	for (i = 0; i < Nelem; i ++) {
		for (n = 0; n < Nca_goal; n ++) cap_goal[i][n] = 0.;
		for (j = 0; j < Ngoal; j ++) {
			n = (int)(goal[j].Celem[i] * (goal[j].Nnebr + 1) + 0.001);
			cap_goal[i][n] ++;
		}
		for (n = 0; n < Nca_goal; n ++) cap_goal[i][n] /= Ngoal;
		printf ("%s\n", element[i]);
		for (n = 0; n < Nca_goal; n ++) printf ("%f\n", cap_goal[i][n]);
	}
	system ("mkdir ../out && mkdir ../out/CE");
	output = fopen ("../out/CE/ce.dat", "w");
	fprintf (output, "C ");
	for (i = 0; i < Nelem; i ++) {
		fprintf (output, "%s ", element[i]);
	}
	fprintf (output, "\n");
	for (n = 0; n < Nca_goal; n ++) {
		fprintf (output, "%.2f ", (double)n / (double)(Nca_goal - 1));
		for (i = 0; i < Nelem; i ++) {
			fprintf (output, "%f ", cap_goal[i][n]);
		}
		fprintf (output, "\n");
	}
	fclose (output);
	for (n = 0; n < Ngoal; n ++) free (goal[n].Celem);
	free (goal);

	printf ("\n\n/********************kMC Cluster mapping*******************/\n\n");
	N = nMol;
	AllocMem (Nelem_atom, Nelem, int);
	for (i = 0; i < Nelem; i ++) Nelem_atom[i] = 0;
	DO_MOL {
		for (i = 0; i < Nelem; i ++) {
			if (strcmp (mol[n].elem, element[i]) == 0) Nelem_atom[i] ++;
		}
	}
	for (i = 1; i < Nelem; i ++) {
		if (Nelem_atom[i] > (int)(nMol * Celem[i])) {
			j = 0;
			while (1) {
				n = rand() % nMol;
				if (strcmp (mol[n].elem, element[i]) == 0) {
					sprintf (mol[n].elem, "%s", element[0]);
					j ++;
				}
				if (j >= Nelem_atom[i] - (int)(nMol * Celem[i])) break;
			}
		}
	}
	for (i = 1; i < Nelem; i ++) {
		if (Nelem_atom[i] < (int)(nMol * Celem[i])) {
			j = 0;
			while (1) {
				n = rand() % nMol;
				if (strcmp (mol[n].elem, element[0]) == 0) {
					sprintf (mol[n].elem, "%s", element[i]);
					j ++;
				}
				if (j >= (int)(nMol * Celem[i]) - Nelem_atom[i]) break;
			}
		}
	}
	free (Nelem_atom);
/*
	DO_MOL {
		sprintf (mol[n].elem, "%s", element[0]);
	}
	for (i = 1; i < Nelem; i ++) {
		j = 0;
		while (1) {
			n = rand() % nMol;
			if (strcmp (mol[n].elem, element[0]) == 0) {
				sprintf (mol[n].elem, "%s", element[i]);
				j ++;
			}
			if (j >= Celem[i] * nMol) break;
		}
	}
	for (i = 0; i < Nelem; i ++) {
		printf ("%s: C = %.3f%%\n", element[i], Celem[i] * 100.);
	}
*/
	//find rmin1, rmin2 and rmin3
	drmin1 = drmin2 = drmin3 = 1.e10;
	for (i = 1; i < nMol; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, mol[j1].r, mol[j2].r);
		PeriodBoundary (Vdr);
		dr = VLen (Vdr);
		if (drmin1 > dr) drmin1 = dr;
	}
	for (i = 1; i < nMol; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, mol[j1].r, mol[j2].r);
		PeriodBoundary (Vdr);
		dr = VLen (Vdr);
		if ((drmin2 > dr) && (dr - drmin1 > 0.01)) drmin2 = dr;
	}
	for (i = 1; i < nMol; i ++) {
		j1 = 0;
		j2 = i;
		VSub (Vdr, mol[j1].r, mol[j2].r);
		PeriodBoundary (Vdr);
		dr = VLen (Vdr);
		if ((drmin3 > dr) && (dr - drmin2 > 0.01)) drmin3 = dr;
	}
	printf ("the nearest proximity distance is %f\nthe secondary proximity distance is %f\n", drmin1, drmin2);
	printf ("the third proximity distance is %f\n", drmin3);

			//build cell list
			if (Ndr == 1) rNebr = drmin1 + 0.01;
			else if (Ndr == 2) rNebr = drmin2 + 0.01;
			else rNebr = drmin3 + 0.01;
			for (i = 0; i < N; i ++) mol[i].Nnebr = 0;
			nebrTabMax = 200 * nMol;
			VSCopy (cells, 1. / rNebr, region);
			printf ("cells: [%d %d %d]\n", cells.x, cells.y, cells.z);
			AllocMem (cellList, VProd (cells) + nMol, int);
			AllocMem (nebrTab, 2 * nebrTabMax, int);
			rrNebr = Sqr (rNebr);
			VDiv (invWid, cells, region);
			for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
			DO_MOL{
//				VSAdd (rs, mol[n].r, 0.5, region);
				VCopy (rs, mol[n].r);
				VMul (cc, rs, invWid);
				c = VLinear (cc, cells) + nMol;
				if (c >= (VProd(cells) + nMol) || c < 0) {
					printf ("error(nebrList.c): cellList is not big enough(measurement.c),\n \
						cc:(%d %d %d), cells(%d %d %d)\n \
						rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
						rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
						region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
					PrintMovie ();
					exit (1);
				}
				cellList[n] = cellList[c];
				cellList[c] = n;
			}
			nebrTabLen = 0;
			for (m1z = 0; m1z < cells.z; m1z ++){
				for (m1y = 0; m1y < cells.y; m1y ++){
					for (m1x = 0; m1x < cells.x; m1x ++){
						VSet (m1v, m1x, m1y, m1z);
						m1 = VLinear (m1v, cells) + nMol;
						for (offset = 0; offset < N_OFFSET; offset ++){
							VAdd (m2v, m1v, vOff[offset]);
							VZero (shift);
							VCellWrapAll ();
							m2 = VLinear (m2v, cells) + nMol;
							DO_CELL (j1, m1){
								DO_CELL (j2,m2){
									if (m1 != m2 || j2 < j1){
										VSub (Vdr, mol[j1].r, mol[j2].r);
										VVSub (Vdr, shift);
										if (VLenSq (Vdr) < rrNebr){
											if (nebrTabLen >= nebrTabMax) {
												printf ("error(measurement.c): ERR_TOO_MANY_NEBRS");
												exit (1);
											}
											nebrTab[2 * nebrTabLen] = j1;
											nebrTab[2 * nebrTabLen + 1] = j2;
											++ nebrTabLen;
											i = mol[j1].Nnebr;
											j = mol[j2].Nnebr;
											mol[j1].nebr[i] = j2;
											mol[j2].nebr[j] = j1;
											mol[j1].Nnebr ++;
											mol[j2].Nnebr ++;
										}
									}
								}
							}
						}
					}
				}
			}

	Nca = mol[0].Nnebr + 2;
	printf ("the atom number of each cluster is %d\n", Nca - 1);
	AllocMem2 (cap, Nelem, Nca, double); //the proportion of each kind of clusters

	if (Nca != Nca_goal) {
		printf ("error: the number of neighbours is different between CONTCAR and goal cluters\n");
		exit (1);
	}
/*
	DO_MOL {
		if (mol[n].Nnebr + 2 != Nca) {
			printf ("error:box region may be wrong\n");
			exit (1);
		}
	}
*/
	loop = 0;
	diff_old = 1.e100;
	while (1) {
		while (1) {
			j1 = rand () % nMol;
			j2 = rand () % nMol;
			if (Norder == 0) {
				if (strcmp(mol[j1].elem, mol[j2].elem) != 0) break;
			} else {
				if (strcmp(mol[j1].elem, mol[j2].elem) != 0 && mol[j1].r.x > region.x * 0.4 && mol[j2].r.x > region.x * 0.4) break;
			}
		}
		strcpy (elem, mol[j1].elem);
		strcpy (mol[j1].elem, mol[j2].elem);
		strcpy (mol[j2].elem, elem);

		for (i = 0; i < Nelem; i ++) {
			for (n = 0; n < Nca; n ++) cap[i][n] = 0.;
			for (j = 0; j < nMol; j ++) {
				n = 0;
				if (strcmp(mol[j].elem, element[i]) == 0) n ++;
				for (k = 0;k < mol[j].Nnebr; k ++) {
					id = mol[j].nebr[k];
					if (strcmp(mol[id].elem, element[i]) == 0) n ++;
				}
				cap[i][n] ++;
			}
			for (n = 0; n < Nca; n ++) cap[i][n] /= (double)nMol;
		}

		diff = 0;
		for (j = 0; j < Nelem; j ++) {
			for (i = 0; i < Nca; i ++) {
				if (error == 1) diff += Sqr (cap[j][i] - cap_goal[j][i]);
				else if (error == 2) diff += fabs (cap[j][i] - cap_goal[j][i]);
			}
		}
		if (error == 1) diff = sqrt(diff / Nca);
		else if (error == 2) diff /= Nca;

		if (diff > diff_old && (double)rand() / RAND_MAX > prob) {
			strcpy (elem, mol[j1].elem);
			strcpy (mol[j1].elem, mol[j2].elem);
			strcpy (mol[j2].elem, elem);
		} else {
			diff_old = diff;
		}

		if (loop % (loop_max / 1000) == 0) printf ("--task: %.1f%%, diff = %e\n", (double)loop / (double)loop_max * 100., diff_old);
		if (loop == loop_max) {
			output = fopen ("../out/CE/kMC_mapping.dat", "a+");
			fprintf (output, "C ");
			for (i = 0; i < Nelem; i ++) {
				fprintf (output, "%s ", element[i]);
				fprintf (output, "%s(goal) ", element[i]);
			}
			fprintf (output, "\n");
			for (i = 0; i < Nca; i ++) {
				fprintf (output, "%.3f ", (double)i / (double)(Nca - 1));
				for (j = 0; j < Nelem; j ++) {
					fprintf (output, "%f %f ",  cap[j][i], cap_goal[j][i]);
				}
				fprintf (output, "\n");
			}
			
			fprintf (output, "--task: %.3f%%, diff = %e\n", (double)loop / (double)loop_max * 100., diff_old);
			fclose (output);
			output = fopen ("../out/CE/kMC_mapping.movie", "a+");
			
			fprintf (output, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3\n", \
				 nMol, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
			
			DO_MOL {
				fprintf (output, "%d %s %f %f %f\n", \
					n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10)); //--A
			}
			fclose (output);
		}

		loop ++;
		if (loop > loop_max) break;
	}
	
	free (element[0]);
	free (element);
	free (cap_goal[0]);
	free (cap_goal);
	free (cap[0]);
	free (cap);
	free (Celem);
}
