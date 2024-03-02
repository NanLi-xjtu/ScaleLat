#include "mc.h"

Atom *atom;
VecI cellSize;
VecR region;
double TUnit, EUnit, VUnit, lUnit, eUnit, T_max, T_min, dT, dT_min, dT_max, threshold, *pairDists, NLOOP, stepMovie, stepAvg;
int *Natoms, Nelems, NeighbourMax, Ndists = 1;
long long Nloop, stepCount;
int nthreads, seed;
char MC_method[128], cell_order[128];
int stepVacancy, stepAtom, movie_count = 0;

void Setup ()
{
	FILE *input, *input1;
	char filename[128], line[1024];
	int i, j, k, n, temp, *count, range;
	Atom atom_exchange;

	//Unit
	TUnit = 300.; //--K
	eUnit = EUnit = TUnit * kB; //--J
	VUnit = EUnit / eleChar; //--eV
	lUnit = 1.e-10; //--m

	//read atom number from super cell
	AllocMem (Natoms, Nelems+1, int);
	sprintf (filename, "../out/superCell/%d_%d_%d/Supercell.xyz", cellSize.x, cellSize.y, cellSize.z);
	input = ReadFile (filename);
	fgets (line, 1024, input);
	sscanf (line, "%d", &Natoms[0]);
	printf ("total atoms: %d\n", Natoms[0]);
	nMol = Natoms[0];
	Nloop = (long long)NLOOP;

	//read region
	sprintf (filename, "../out/superCell/%d_%d_%d/Lattice_supercell", cellSize.x, cellSize.y, cellSize.z);
	input1 = ReadFile (filename);
	fgets (line, 1024, input1);
	sscanf (line, "%lg", &region.x);
	fgets (line, 1024, input1);
	sscanf (line, "%lg", &region.y);
	fgets (line, 1024, input1);
	sscanf (line, "%lg", &region.z);
	fclose (input1);

	//allocate memory
	AllocMem (atom, Natoms[0], Atom);
	DO_ATOM {
		AllocMem (atom[n].nebrLists, Ndists, NebrLists);
		for (i = 0; i < Ndists; i ++) AllocMem (atom[n].nebrLists[i].nebr, NeighbourMax, int);
	}
	AllocMem (pairDists, Ndists, double);

	//init coordination and element
	fgets (line, 1024, input);
	n = 0;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;
		sscanf (line, "%s %lg %lg %lg\n", atom[n].element, &atom[n].r.x, &atom[n].r.y, &atom[n].r.z);
		n ++;
	}
	DO_ATOM {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp(atom[n].element, elemName[i]) == 0) atom[n].type = i+1;
		}
	}
	fclose (input);

	//atom numbers of different elements
	printf ("atom number: ");
	for (i = 1; i <= Nelems; i ++) {
		Natoms[i] = 0;
		DO_ATOM {
			if (strcmp(atom[n].element, elemName[i-1]) == 0) Natoms[i] ++;
		}
		printf ("%s %d %.2f%% ", elemName[i-1], Natoms[i], (double)(Natoms[i])/(double)(Natoms[0])*100.);
	}
	printf ("\n");

	//the initial elements are randomly distrubuted
	if (strcmp(cell_order, "disorder") == 0) {
		printf ("\natom strcuture is disorder\n");
		DO_ATOM {
			i = rand () % Natoms[0];
			atom_exchange = atom[n];
			atom[n] = atom[i];
			atom[i] = atom_exchange;
		}
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		DO_ATOM {
			j = rand () % Natoms[0];
			range = 0;
			for (i = 0; i < Nelems; i ++) {
				if (j >= range && j < range + Natoms[i+1]) break;
				range += Natoms[i+1];
			}

			if (count[i+1] >= Natoms[i+1]) {
				n --;
				continue;
			}
			strcpy (atom[n].element, elemName[i]);
			count[i+1] ++;
			count[0] ++;
		}
		free (count);
	} else if (strcmp(cell_order, "order") == 0) {
		printf ("\natom strcuture is order\n");
	} else {
		printf ("error(mc.c): cell_order format\n");
		exit (1);
	}
}

void SetupMD ()
{
	int n, i, j, k;

	if (strcmp (mode, "kMC") == 0) {
		printf ("\n-- init eam potential...\n");
		coeff_alloy(filename_EAM);
		init_style_alloy ();
		rCut = cutmax; //--A
	} else rCut = 6.0; //--A

	rCut /= lUnit * 1.e10; //--dimensionless
	rNebrShell = 0.; //--A
	rNebrShell /= lUnit * 1.e10; //--dimensionless
	printf ("rCut: %f A\n", rCut);
	VScale (region, 1. / (lUnit * 1.e10)); //--dimensionless
	nebrTabMax = nebrTabFac * nMol;
	AllocMem (nebrTab, 2 * nebrTabMax, int);
	AllocMem (mol, nMol, Mol);
	DO_MOL AllocMem (mol[n].nebr_id, nebrTabFac, int);
	DO_MOL {
		VSCopy (mol[n].r, 1. / (lUnit * 1.e10), atom[n].r); //dismensionless
		sprintf (mol[n].elem, "%s", atom[n].element);
	}

	printf ("\n\n-- build neighbour list...\n\n");
	VSCopy (cells, 1. / (rCut + rNebrShell), region);
	AllocMem (cellList, VProd (cells) + nMol, int);
	BuildNebrList ();
	free (cellList);
//	CheckNebrList ();
	VSCopy (cells, 1. / (rCut), region);
	AllocMem (cellList, VProd (cells) + nMol, int);
	GetNebrs (rCut);
	free (cellList);

	//get atom pair distance
	Get_pair_distance ();
	//get i-th sphere neighbour atoms and pairs lists
	Get_neighbour ();
//	PrintNeighbour ();

	//RDF
	rangeRdf /= lUnit * 1.e10; //dimensionless
	AllocMem2 (histRdf, (int)(Nelems*Nelems+1.), sizeHistRdf, double);
	countRdf = 0;
}

void MonteCarloVacancy (double T)
{
	long long i, elemCount[10]={0};
	int n, vacancy, j, id, nebr, min, count, move_flag = 0, m1, m2;
	double dE, dE2, E_mi, dE_mi, mu, Gamma_AV, dt, t_tot;
	double start, end;
	FILE *output;
	char filename[128], element[10];

	vacancy = -1;
	DO_MOL {
		if (strcmp (mol[n].elem, "vacancy") == 0) vacancy = n;
	}
	if (vacancy == -1) {
		vacancy = rand() % Natoms[0]; //pick one element as a electron vacancy
		strcpy (mol[vacancy].elem, "vacancy");
	}
	count = 0;
	t_tot = 0.;
	start = omp_get_wtime();
	T /= TUnit; //--dimensionless
	for (i = 0; i < Nloop; i ++) {
		stepCount = i;
		if (stepCount % stepRdf == 0) EvalRdf (); //RDF
		if (i % (long long)(stepMovie) == 0) PrintMovie ();
		if (i % (int)stepAvg == 0) {
			//output move steps and real time
//			for (j = 0; j < Nelems; j ++) printf ("%s %.3e ", elemName[j], (double)elemCount[j]);
			printf ("t: %.3e s, ", t_tot);
			sprintf (filename, "../out/real_time.dat");
			output = WriteFile (filename);
			if (move_flag == 0) {
				fprintf (output, "step ");
				for (j = 0; j < Nelems; j ++) fprintf (output, "%s ", elemName[j]);
				fprintf (output, "time(s)\n");
				move_flag = 1;
			}
			fprintf (output, "%e ", (double)i);
			for (j = 0; j < Nelems; j ++) fprintf (output, "%e ", (double)elemCount[j]);
			fprintf (output, "%e\n", t_tot);
			fclose (output);
			//update atom sort
			DO_ATOM {
				sprintf (atom[n].element, "%s", mol[n].elem);
				for (j = 0; j < Nelems; j ++) {
					if (strcmp(atom[n].element, elemName[j]) == 0) atom[n].type = j+1;
				}
			}
			sprintf (filename, "../out/SRO.dat");
			ShortRangeOrder (elemName[0], elemName[1], filename);
			//compute system energy
			ComputeForcesEamPoten_alloy ();
			printf ("<u>: %.4f eV, ", (uSum / nMol) * VUnit);
			end = omp_get_wtime();
			printf ("step: %.2e(%.4f%%), time: %5.2f s\n", (double)i, (double)i / (double)Nloop * 100., (end - start));
			start = omp_get_wtime();
		}

		id = rand () % (atom[vacancy].nebrLists[0].Nnebrs); //pick one neighbour to replace the vacancy
		nebr = atom[vacancy].nebrLists[0].nebr[id];
		dE = Get_dE_new (vacancy, nebr); //--eV
//		dE2 = Get_dE (vacancy, nebr); //--eV
//		printf ("%f %f eV\n", dE, dE2);getchar ();

		dE /= (eUnit / eleChar); //--dimensionless
//		if (count > stepVacancy) printf ("there is a random exchange in case endless loop\n");
		if (dE < 0 || (double)rand() / (double)RAND_MAX < exp(- dE / T) || count > stepVacancy) {
			//calculate real time
			E_mi = 0.;
			for (j = 0; j < Nelems; j ++) {
				if (strcmp (mol[nebr].elem, elemName[j]) == 0) elemCount[j] ++;
				if (strcmp (mol[nebr].elem, "Fe") == 0) E_mi = 0.64; //--eV
				if (strcmp (mol[nebr].elem, "Cr") == 0) E_mi = 0.57; //--eV
			}
			E_mi /= (eUnit / eleChar); //--dimensionless
			dE_mi = E_mi - 0.5 * dE;
			mu = 1.e12; //--HZ
			Gamma_AV = mu * exp(- dE_mi / T); //--HZ
			dt = 1. / Gamma_AV * nMol; //--s
			t_tot += dt; //--s

			strcpy (mol[vacancy].elem, mol[nebr].elem);
			strcpy (mol[nebr].elem, "vacancy");
			vacancy = nebr;
			count = 0;
		} else strcpy (mol[vacancy].elem, "vacancy");

		count ++;
	}
}

void MonteCarloExchangeAtoms (double T)
{
	long long i, elemCount = 0;
	int n, vacancy, j, id, nebr, min, count, move_flag = 0, m1, m2;
	double dE;
	double start, end;
	FILE *output;
	char filename[128], element[10];

	start = omp_get_wtime();
	T /= TUnit; //--dimensionless
	for (i = 0; i < Nloop; i ++) {
		stepCount = i;
		if (stepCount % stepRdf == 0) EvalRdf (); //RDF
		if (i % (long long)(stepMovie) == 0) PrintMovie ();
		if (i % (int)stepAvg == 0) {
			//count the atoms exchange step
			printf ("exchange_count: %.3e ", (double)elemCount);
			sprintf (filename, "../out/move_count.dat");
			output = WriteFile (filename);
			if (move_flag == 0) {
				fprintf (output, "step exchange_count\n");
				move_flag = 1;
			}
			fprintf (output, "%e %e\n", (double)i, (double)elemCount);
			fclose (output);
			//update atom sort
			DO_ATOM {
				sprintf (atom[n].element, "%s", mol[n].elem);
				for (j = 0; j < Nelems; j ++) {
					if (strcmp(atom[n].element, elemName[j]) == 0) atom[n].type = j+1;
				}
			}
			sprintf (filename, "../out/SRO.dat");
			ShortRangeOrder (elemName[0], elemName[1], filename);
			//compute system energy
			ComputeForcesEamPoten_alloy ();
			printf ("<u>: %.4f eV ", (uSum / nMol) * VUnit);
			end = omp_get_wtime();
			printf ("%.2e %.4f%% time: %.2f s\n", (double)i, (double)i / (double)Nloop * 100., (end - start));
			start = omp_get_wtime();
		}

		while (1) {
			m1 = rand () % nMol;
			m2 = rand () % nMol;
			if (strcmp(mol[m1].elem, mol[m2].elem) != 0) break;
		}
		dE = Get_dE_AtomsExchange (m1, m2); //--eV, the element of m1 and m2 has been exchanged

		dE /= (eUnit / eleChar); //--dimensionless
		if (dE < 0 || (double)rand() / (double)RAND_MAX < exp(- dE / T)) {
			elemCount ++;
		} else {
			strcpy (element, mol[m1].elem);
			strcpy (mol[m1].elem, mol[m2].elem);
			strcpy (mol[m2].elem, element);
		}
	}
}

void MonteCarloMixed (double T)
{
	long long i, elemCount[10]={0}, exchangeCount;
	int n, vacancy, j, id, nebr, min, count, move_flag = 0, m1, m2;
	double dE;
	double start, end;
	FILE *output;
	char filename[128], element[10];

	vacancy = 0;
	DO_MOL {
		if (strcmp (mol[n].elem, "vacancy") == 0) vacancy = n;
	}
	if (vacancy == 0) {
		vacancy = rand() % Natoms[0]; //pick one element as a electron vacancy
		strcpy (mol[vacancy].elem, "vacancy");
	}
	count = 0;
	start = omp_get_wtime();
	exchangeCount = 0;
	T /= TUnit; //--dimensionless
	for (i = 0; i < Nloop; i ++) {
		stepCount = i;
		if (stepCount % stepRdf == 0) EvalRdf (); //RDF
		if (stepCount % (long long)(stepMovie) == 0) PrintMovie ();
		if (stepCount % (int)stepAvg == 0) {
			for (j = 0; j < Nelems; j ++) printf ("%s %.3e ", elemName[j], (double)elemCount[j]);
			printf ("atomexchange %.3e ", (double)exchangeCount);
			sprintf (filename, "../out/move_count.dat");
			output = WriteFile (filename);
			if (move_flag == 0) {
				fprintf (output, "step ");
				for (j = 0; j < Nelems; j ++) fprintf (output, "%s ", elemName[j]);
				fprintf (output, "exchangeCount\n");
				move_flag = 1;
			}
			fprintf (output, "%e ", (double)i);
			for (j = 0; j < Nelems; j ++) fprintf (output, "%e ", (double)elemCount[j]);
			fprintf (output, "%e\n", (double)exchangeCount);
			fclose (output);
			//update atom sort
			DO_ATOM {
				sprintf (atom[n].element, "%s", mol[n].elem);
				for (j = 0; j < Nelems; j ++) {
					if (strcmp(atom[n].element, elemName[j]) == 0) atom[n].type = j+1;
				}
			}
			sprintf (filename, "../out/SRO.dat");
			ShortRangeOrder (elemName[0], elemName[1], filename);
			//compute system energy
			ComputeForcesEamPoten_alloy ();
			printf ("<u>: %.4f eV ", (uSum / nMol) * VUnit);
			end = omp_get_wtime();
			printf ("step %.2e (%.4f%%) time: %5.2f s\n", (double)i, (double)i / (double)Nloop * 100., (end - start));
			start = omp_get_wtime();
		}

		id = rand () % (atom[vacancy].nebrLists[0].Nnebrs); //pick one neighbour to replace the vacancy
		nebr = atom[vacancy].nebrLists[0].nebr[id];
		dE = Get_dE_new (vacancy, nebr); //--eV

		dE /= (eUnit / eleChar); //--dimensionless
		if (dE < 0 || (double)rand() / (double)RAND_MAX < exp(- dE / T)) {
			for (j = 0; j < Nelems; j ++) {
				if (strcmp (mol[nebr].elem, elemName[j]) == 0) elemCount[j] ++;
			}
			strcpy (mol[vacancy].elem, mol[nebr].elem);
			strcpy (mol[nebr].elem, "vacancy");
			vacancy = nebr;
			count = 0;
		} else if (count > stepVacancy) {       //if the vancancy moving get stuck for more than a certain steps, 
			for (j = 0; j < Nelems; j ++) { //exchanging atoms will be performed for a certain steps.
				if (strcmp (mol[nebr].elem, elemName[j]) == 0) elemCount[j] ++;
			}
			strcpy (mol[vacancy].elem, mol[nebr].elem);
			strcpy (mol[nebr].elem, "vacancy");
			vacancy = nebr;
			for (j = 0; j < stepAtom; j ++) {
				while (1) {
					m1 = rand () % nMol;
					m2 = rand () % nMol;
					if (strcmp(mol[m1].elem, mol[m2].elem) != 0 && strcmp(mol[m1].elem, "vacancy") != 0 && \
					    strcmp(mol[m2].elem, "vacancy") != 0) break;
				}
				dE = Get_dE_AtomsExchange (m1, m2); //--eV, the element of m1 and m2 has been exchanged

				dE /= (eUnit / eleChar); //--dimensionless
				if (dE < 0 || (double)rand() / (double)RAND_MAX < exp(- dE / T)) {
					exchangeCount ++;
				} else {
					strcpy (element, mol[m1].elem);
					strcpy (mol[m1].elem, mol[m2].elem);
					strcpy (mol[m2].elem, element);
				}
			}
			count = 0;
		} else strcpy (mol[vacancy].elem, "vacancy");

		count ++;
	}
}

void Restart ()
{
	FILE *input, *input1;
	char filename[128], line[1024];
	int i, j, k, n, *count, range;
	double temp;

	//Unit
	TUnit = 300.; //--K
	EUnit = TUnit * kB; //--J
	VUnit = EUnit / eleChar; //--eV
	lUnit = 1.e-10; //--m
	eUnit = EUnit; //--J

	//read restart.xyz
	AllocMem (Natoms, Nelems+1, int);
	sprintf (filename, "../in/restart.xyz");
	input = ReadFile (filename);
	fgets (line, 1024, input);
	sscanf (line, "%d", &Natoms[0]);
	printf ("total atoms: %d\n", Natoms[0]);
	nMol = Natoms[0];
	Nloop = (long long)NLOOP;
	fgets (line, 1024, input);
	sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", \
		&region.x, &temp, &temp, &temp, &region.y, &temp, &temp, &temp, &region.z);

	//allocate memory
	AllocMem (atom, Natoms[0], Atom);
	DO_ATOM {
		AllocMem (atom[n].nebrLists, Ndists, NebrLists);
		for (i = 0; i < Ndists; i ++) AllocMem (atom[n].nebrLists[i].nebr, NeighbourMax, int);
	}
	AllocMem (pairDists, Ndists, double);

	//read coordination and element
	n = 0;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;
		sscanf (line, "%lg %s %lg %lg %lg\n", &temp, atom[n].element, &atom[n].r.x, &atom[n].r.y, &atom[n].r.z);
		n ++;
	}
	DO_ATOM {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp(atom[n].element, elemName[i]) == 0) atom[n].type = i+1;
		}
	}
	fclose (input);

	//atom numbers of different elements
	printf ("atom number: ");
	for (i = 1; i <= Nelems; i ++) {
		Natoms[i] = 0;
		DO_ATOM {
			if (strcmp(atom[n].element, elemName[i-1]) == 0) Natoms[i] ++;
		}
		printf ("%s %d %.2f%% ", elemName[i-1], Natoms[i], (double)(Natoms[i])/(double)(Natoms[0])*100.);
	}
	printf ("\n");
}

double Get_dE (int vacancy, int nebr)
{
	int i;
	double uSum0, uSum1;
	VecR vacancyR;

	ComputeForcesEamPoten_alloy ();
	uSum0 = uSum;

	//move the vacancy
	strcpy (mol[vacancy].elem, mol[nebr].elem);
	strcpy (mol[nebr].elem, "vacancy");
	ComputeForcesEamPoten_alloy ();
	uSum1 = uSum;

	strcpy (mol[nebr].elem, mol[vacancy].elem);
	strcpy (mol[vacancy].elem, "vacancy");

	return ((uSum1 - uSum0) * VUnit); //--eV
}

double Get_dE_vacancy (int vacancy)
{
	int n1, n2, i, j, k, itype, jtype, m, n;
	double F, F1, Fsum, F1sum, dphi, dE;
	VecR drVec, DR_A;
	double r, rsq, p, *coeff, Fp, z2, recip, phi, rho;

	//calculate embedding energy F
	Fsum = 0.;
	#pragma omp parallel for reduction(+:Fsum) private(i,itype,k,rho,n2,j,drVec,DR_A,rsq,jtype,p,m,coeff,Fp,F) num_threads(nthreads)
	for (n1 = 0; n1 < mol[vacancy].Nnebr + 1; n1 ++) {
		if (n1 < mol[vacancy].Nnebr) i = mol[vacancy].nebr_id[n1];
		else i = vacancy;
		itype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
		}
		if (itype == 0) {
			printf ("error1(mc.c):itype\n");
			exit (1);
		}
		rho = 0.;
		for (n2 = 0; n2 < mol[i].Nnebr; n2 ++) {
			j = mol[i].nebr_id[n2];
			VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
			VWrapAll (drVec);
			VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
			rsq = VLenSq (DR_A); //--A2
			jtype = 0;
			for (k = 1; k <= ntypes; k ++) {
				if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
			}
			if (jtype == 0) {
				printf ("error2(mc.c):jtype\n");
				exit (1);
			}
			p = sqrt(rsq)*rdr + 1.0;
			m = (int) (p);
			m = MIN(m,nr-1);
			p -= m;
			p = MIN(p,1.0);
			coeff = rhor_spline[type2rhor[jtype][itype]][m];
      			rho += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		}
		p = rho*rdrho + 1.0;
		m = (int) (p);
		m = MAX(1,MIN(m,nrho-1));
		p -= m;
		p = MIN(p,1.0);
		coeff = frho_spline[type2frho[itype]][m];
		Fp = (coeff[0]*p + coeff[1])*p + coeff[2];
		F = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		if (rho > rhomax) F += Fp * (rho-rhomax);
		Fsum += F; //--eV
	}
	//consider the vacancy, calculate embedding energy F1
	F1sum = 0.;
	#pragma omp parallel for reduction(+:F1sum) private(i,itype,k,rho,n2,j,drVec,DR_A,rsq,jtype,p,m,coeff,Fp,F1) num_threads(nthreads)
	for (n1 = 0; n1 < mol[vacancy].Nnebr; n1 ++) {
		i = mol[vacancy].nebr_id[n1];
		itype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
		}
		if (itype == 0) {
			printf ("error3(mc.c):itype\n");
			exit (1);
		}
		rho = 0.;
		for (n2 = 0; n2 < mol[i].Nnebr; n2 ++) {
			j = mol[i].nebr_id[n2];
			if (j == vacancy) continue;
			VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
			VWrapAll (drVec);
			VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
			rsq = VLenSq (DR_A); //--A2
			jtype = 0;
			for (k = 1; k <= ntypes; k ++) {
				if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
			}
			if (jtype == 0) {
				printf ("error4(mc.c):jtype\n");
				exit (1);
			}
			p = sqrt(rsq)*rdr + 1.0;
			m = (int) (p);
			m = MIN(m,nr-1);
			p -= m;
			p = MIN(p,1.0);
			coeff = rhor_spline[type2rhor[jtype][itype]][m];
      			rho += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		}
		p = rho*rdrho + 1.0;
		m = (int) (p);
		m = MAX(1,MIN(m,nrho-1));
		p -= m;
		p = MIN(p,1.0);
		coeff = frho_spline[type2frho[itype]][m];
		Fp = (coeff[0]*p + coeff[1])*p + coeff[2];
		F1 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		if (rho > rhomax) F1 += Fp * (rho-rhomax);
		F1sum += F1; //--eV
	}
	//consider the vacancy, calculate pire potential
	dphi = 0.;
	i = vacancy;
	for (n1 = 0; n1 < mol[vacancy].Nnebr; n1 ++) {
		j = mol[vacancy].nebr_id[n1];
		itype = jtype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
			if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
		}
		if (itype == 0 || jtype == 0) {
			printf ("error5(mc.c):itype jtype\n");
			exit (1);
		}
		VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
		VWrapAll (drVec);
		VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
		rsq = VLenSq (DR_A); //--A2
		r = sqrt(rsq);
		p = r*rdr + 1.0;
		m = (int) (p);
		m = MIN(m,nr-1);
		p -= m;
		p = MIN(p,1.0);
		coeff = z2r_spline[type2z2r[itype][jtype]][m];
		z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		recip = 1.0/r;
		phi = z2*recip; //--eV
		dphi += phi; //--eV
	}
	dE = Fsum - F1sum + dphi; //--eV

	return (dE);
}

double Get_dE_new (int vacancy, int nebr)
{
	double dE1, dE2, dE;

	strcpy (mol[vacancy].elem, mol[nebr].elem);
	dE1 = Get_dE_vacancy (vacancy);
	dE2 = Get_dE_vacancy (nebr);
	dE = dE1 - dE2;
	strcpy (mol[vacancy].elem, "vacancy");

	return (dE);
}

double Get_dE_OneAtomChange (int id, char *element)
{
	int n1, n2, i, j, k, itype, jtype, m, n;
	double F, F1, Fsum, F1sum, dphi, dphi1, dE;
	VecR drVec, DR_A;
	double r, rsq, p, *coeff, Fp, z2, recip, phi, rho;

	//calculate embedding energy F
	Fsum = 0.;
	#pragma omp parallel for reduction(+:Fsum) private(i,itype,k,rho,n2,j,drVec,DR_A,rsq,jtype,p,m,coeff,Fp,F) num_threads(nthreads)
	for (n1 = 0; n1 < mol[id].Nnebr + 1; n1 ++) {
		if (n1 < mol[id].Nnebr) i = mol[id].nebr_id[n1];
		else i = id;
		itype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
		}
		if (strcmp(mol[i].elem, "vacancy") == 0) continue;
		if (itype == 0) {
			printf ("error6(mc.c):itype\n");
			exit (1);
		}
		rho = 0.;
		for (n2 = 0; n2 < mol[i].Nnebr; n2 ++) {
			j = mol[i].nebr_id[n2];
			VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
			VWrapAll (drVec);
			VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
			rsq = VLenSq (DR_A); //--A2
			jtype = 0;
			for (k = 1; k <= ntypes; k ++) {
				if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
			}
			if (strcmp(mol[j].elem, "vacancy") == 0) continue;
			if (jtype == 0) {
				printf ("error7(mc.c):jtype %s\n", mol[j].elem);
				exit (1);
			}
			p = sqrt(rsq)*rdr + 1.0;
			m = (int) (p);
			m = MIN(m,nr-1);
			p -= m;
			p = MIN(p,1.0);
			coeff = rhor_spline[type2rhor[jtype][itype]][m];
      			rho += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		}
		p = rho*rdrho + 1.0;
		m = (int) (p);
		m = MAX(1,MIN(m,nrho-1));
		p -= m;
		p = MIN(p,1.0);
		coeff = frho_spline[type2frho[itype]][m];
		Fp = (coeff[0]*p + coeff[1])*p + coeff[2];
		F = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		if (rho > rhomax) F += Fp * (rho-rhomax);
		Fsum += F; //--eV
	}
	//calculate pire potential
	dphi = 0.;
	i = id;
	for (n1 = 0; n1 < mol[id].Nnebr; n1 ++) {
		j = mol[id].nebr_id[n1];
		itype = jtype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
			if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
		}
		if (strcmp(mol[i].elem, "vacancy") == 0 || strcmp(mol[j].elem, "vacancy") == 0) continue;
		if (itype == 0 || jtype == 0) {
			printf ("error8(mc.c):itype jtype\n");
			exit (1);
		}
		VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
		VWrapAll (drVec);
		VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
		rsq = VLenSq (DR_A); //--A2
		r = sqrt(rsq);
		p = r*rdr + 1.0;
		m = (int) (p);
		m = MIN(m,nr-1);
		p -= m;
		p = MIN(p,1.0);
		coeff = z2r_spline[type2z2r[itype][jtype]][m];
		z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		recip = 1.0/r;
		phi = z2*recip; //--eV
		dphi += phi; //--eV
	}

	//change the element of atom[id]
	strcpy (mol[id].elem, element);

	//when id is changed to element, calculate embedding energy F1
	F1sum = 0.;
	#pragma omp parallel for reduction(+:F1sum) private(i,itype,k,rho,n2,j,drVec,DR_A,rsq,jtype,p,m,coeff,Fp,F) num_threads(nthreads)
	for (n1 = 0; n1 < mol[id].Nnebr + 1; n1 ++) {
		if (n1 < mol[id].Nnebr) i = mol[id].nebr_id[n1];
		else i = id;
		itype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
		}
		if (strcmp(mol[i].elem, "vacancy") == 0) continue;
		if (itype == 0) {
			printf ("error9(mc.c):itype\n");
			exit (1);
		}
		rho = 0.;
		for (n2 = 0; n2 < mol[i].Nnebr; n2 ++) {
			j = mol[i].nebr_id[n2];
			VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
			VWrapAll (drVec);
			VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
			rsq = VLenSq (DR_A); //--A2
			jtype = 0;
			for (k = 1; k <= ntypes; k ++) {
				if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
			}
			if (strcmp(mol[j].elem, "vacancy") == 0) continue;
			if (jtype == 0) {
				printf ("error10(mc.c):jtype\n");
				exit (1);
			}
			p = sqrt(rsq)*rdr + 1.0;
			m = (int) (p);
			m = MIN(m,nr-1);
			p -= m;
			p = MIN(p,1.0);
			coeff = rhor_spline[type2rhor[jtype][itype]][m];
      			rho += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		}
		p = rho*rdrho + 1.0;
		m = (int) (p);
		m = MAX(1,MIN(m,nrho-1));
		p -= m;
		p = MIN(p,1.0);
		coeff = frho_spline[type2frho[itype]][m];
		Fp = (coeff[0]*p + coeff[1])*p + coeff[2];
		F = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		if (rho > rhomax) F += Fp * (rho-rhomax);
		F1sum += F; //--eV
	}
	//when id is changed to element, calculate pire potential
	dphi1 = 0.;
	i = id;
	for (n1 = 0; n1 < mol[id].Nnebr; n1 ++) {
		j = mol[id].nebr_id[n1];
		itype = jtype = 0;
		for (k = 1; k <= ntypes; k ++) {
			if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
			if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
		}
		if (strcmp(mol[i].elem, "vacancy") == 0 || strcmp(mol[j].elem, "vacancy") == 0) continue;
		if (itype == 0 || jtype == 0) {
			printf ("error11(mc.c):itype jtype\n");
			exit (1);
		}
		VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
		VWrapAll (drVec);
		VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
		rsq = VLenSq (DR_A); //--A2
		r = sqrt(rsq);
		p = r*rdr + 1.0;
		m = (int) (p);
		m = MIN(m,nr-1);
		p -= m;
		p = MIN(p,1.0);
		coeff = z2r_spline[type2z2r[itype][jtype]][m];
		z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		recip = 1.0/r;
		phi = z2*recip; //--eV
		dphi1 += phi; //--eV
	}
	dE = F1sum - Fsum + dphi1 - dphi; //--eV

	return (dE);
}

double Get_dE_AtomsExchange (int m1, int m2)
{
	double dE1, dE2, dE;
	char element[10];

	strcpy (element, mol[m1].elem);
	dE1 = Get_dE_OneAtomChange (m1, mol[m2].elem);
	dE2 = Get_dE_OneAtomChange (m2, element);
	dE = dE1 + dE2;

	return (dE); //--eV
}

void Get_pair_distance ()
{
	int m, n, i, j, k;
	VecR dr;
	double distMin, threshold1;
	FILE *input;
	char line[1024];

	printf ("\n-- getting atom pair distance...\n");

	//get the nearest atomic distances
	for (i = 0; i < Ndists; i ++) {
		distMin = 1.e100;
		if (i == 0) threshold1 = 0.;
		else threshold1 = pairDists[i-1];
		for (k = 0; k < nebrTabLen; k ++) {
			m = nebrTab[2 * k];
			n = nebrTab[2 * k + 1];
			dr.x = fabs(atom[m].r.x - atom[n].r.x);
			dr.y = fabs(atom[m].r.y - atom[n].r.y);
			dr.z = fabs(atom[m].r.z - atom[n].r.z);
			PeriodBoundary (dr);
			if (VLen (dr) < distMin && fabs(VLen (dr)) > threshold + threshold1) {
				distMin = fabs(VLen (dr));
			}
		}
		printf ("The %d_nearest atomic distance is: %e A\n", i+1, distMin);
		pairDists[i] = distMin; //--A
		if (pairDists[i] > (rCut+rNebrShell)*lUnit*1.e10) {
			printf ("error(mc.c): rCut is too small\n");
			exit (1);
		}
	}
}

void Get_neighbour ()
{
	int i, n, j, m, k, j1, j2;
	VecR dr;

	printf ("\n-- getting atom neighbours...\n");
	//get the neighbour atoms
	DO_ATOM {
		for (i = 0; i < Ndists; i ++) {
			atom[n].nebrLists[i].Nnebrs = 0;
		}
	}
	for (i = 0; i < Ndists; i ++) {
		for (k = 0; k < nebrTabLen; k ++) {
			m = nebrTab[2 * k];
			n = nebrTab[2 * k + 1];
			dr.x = fabs(atom[m].r.x - atom[n].r.x);
			dr.y = fabs(atom[m].r.y - atom[n].r.y);
			dr.z = fabs(atom[m].r.z - atom[n].r.z);
			PeriodBoundary (dr);
			if (fabs(VLen (dr) - pairDists[i]) < threshold) {
				j1 = atom[m].nebrLists[i].Nnebrs;
				j2 = atom[n].nebrLists[i].Nnebrs;
				atom[m].nebrLists[i].nebr[j1] = n;
				atom[m].nebrLists[i].Nnebrs ++;
				atom[n].nebrLists[i].nebr[j2] = m;
				atom[n].nebrLists[i].Nnebrs ++;
			}
		}
		if (pairDists[i] > (rCut+rNebrShell)*lUnit*1.e10) {
			printf ("error(mc.c): rCut is too small\n");
			exit (1);
		}
	}
	printf ("ok\n");
}

void MonteCarloDoubleT (double temperature)
{
	long long i, elemCount = 0;
	int n, vacancy, j, id, nebr, min, count, move_flag = 0, m1, m2;
	double dE, sro, sro_old, T;
	double start, end;
	FILE *output;
	char filename[128], element[10];

	start = omp_get_wtime();
	T = 200.; //--K
	T /= TUnit; //--dimensionless
	for (i = 0; i < Nloop; i ++) {
		stepCount = i;
		if (stepCount % stepRdf == 0) EvalRdf (); //RDF
		if (i % (long long)(stepMovie) == 0) PrintMovie ();
		if (i % 1000 == 0) {
			//count the atoms exchange step
			printf ("exchange_count: %.3e ", (double)elemCount);
			sprintf (filename, "../out/move_count.dat");
			output = WriteFile (filename);
			if (move_flag == 0) {
				fprintf (output, "step exchange_count\n");
				move_flag = 1;
			}
			fprintf (output, "%e %e\n", (double)i, (double)elemCount);
			fclose (output);
			//update atom sort
			DO_ATOM {
				sprintf (atom[n].element, "%s", mol[n].elem);
				for (j = 0; j < Nelems; j ++) {
					if (strcmp(atom[n].element, elemName[j]) == 0) atom[n].type = j+1;
				}
			}
			sprintf (filename, "../out/SRO.dat");
			sro = ShortRangeOrder (elemName[0], elemName[1], filename);
			if (stepCount == 0) sro_old = sro;
			//compute system energy
			ComputeForcesEamPoten_alloy ();
			printf ("<u>: %.4f eV ", (uSum / nMol) * VUnit);
			end = omp_get_wtime();
			printf ("%.2e %.4f%% time: %.2f s\n", (double)i, (double)i / (double)Nloop * 100., (end - start));
			if (fabs(sro - sro_old) > 0.05) break;
			start = omp_get_wtime();
		}

		while (1) {
			m1 = rand () % nMol;
			m2 = rand () % nMol;
			if (strcmp(mol[m1].elem, mol[m2].elem) != 0 && strcmp(mol[m1].elem, "vacancy") != 0 && \
			    strcmp(mol[m2].elem, "vacancy") != 0) break;
		}
		dE = Get_dE_AtomsExchange (m1, m2); //--eV, the element of m1 and m2 has been exchanged

		dE /= (eUnit / eleChar); //--dimensionless
		if (dE < 0 || (double)rand() / (double)RAND_MAX < exp(- dE / T)) {
			elemCount ++;
		} else {
			strcpy (element, mol[m1].elem);
			strcpy (mol[m1].elem, mol[m2].elem);
			strcpy (mol[m2].elem, element);
		}
	}
	MonteCarloVacancy (temperature);
}

void PrintAtoms ()
{
	int n;
	FILE *output;
	char filename[128];

	sprintf (filename, "../out/atoms.movie");
	output = WriteFile (filename);
	fprintf (output, "%d\nSolutionReader properties=id:I:1:pos:R:3:type:I:1\n", Natoms[0]);
	DO_ATOM fprintf (output, "%6d%20.10f%20.10f%20.10f %d %s\n", n, atom[n].r.x, atom[n].r.y, atom[n].r.z, atom[n].type, atom[n].element);

	fclose (output);
}

void PrintMovie ()
{
	char filename[128];
	FILE *movie, *restart, *xyz;
	int n;

	if ((movie = fopen ("../out/KMC/out.movie", "a+")) == NULL){
		printf ("\nopen movie file error");
		getchar ();
		exit (1);
	}
	if ((restart = fopen ("../out/restart.xyz", "w")) == NULL){
		printf ("\nopen restart file error");
		getchar ();
		exit (1);
	}
	sprintf (filename, "../out/KMC/out_%d.xyz", movie_count);
	if ((xyz = fopen (filename, "a+")) == NULL){
		printf ("\nopen xyz file error");
		getchar ();
		exit (1);
	}

	fprintf (movie, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3\n", \
		 nMol, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
	DO_MOL {
		fprintf (movie, "%d %s %f %f %f\n", \
			n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10)); //--A
	}

	fprintf (xyz, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3\n", \
		 nMol, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
	DO_MOL {
		fprintf (xyz, "%d %s %f %f %f\n", \
			n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10)); //--A
	}

	fprintf (restart, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3\n", \
		 nMol, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));
	DO_MOL {
		fprintf (restart, "%d %s %f %f %f\n", \
			n, mol[n].elem, mol[n].r.x * (lUnit * 1.e10), mol[n].r.y * (lUnit * 1.e10), mol[n].r.z * (lUnit * 1.e10)); //--A
	}

	movie_count ++;
	fclose (movie);
	fclose (restart);
	fclose (xyz);
}

void PrintNeighbour ()
{
	int m, n, i, nebr;
	FILE *fp;
	char filename[128];

	sprintf (filename, "../out/neighbers.dat");
	fp = WriteFile (filename);
	DO_ATOM {
		fprintf (fp, "%5d%20.10f%20.10f%20.10f %s\n", n, atom[n].r.x, atom[n].r.y, atom[n].r.z, atom[n].element);
		for (i = 0; i < Ndists; i ++) {
			fprintf (fp, "distance = %e:\n", pairDists[i]);
			for (m = 0; m < atom[n].nebrLists[i].Nnebrs; m ++) {
				nebr = atom[n].nebrLists[i].nebr[m];
				fprintf (fp, "%5d%20.10f%20.10f%20.10f %s\n", nebr, atom[nebr].r.x, atom[nebr].r.y, atom[nebr].r.z, atom[nebr].element);
			}
		}
		fprintf (fp, "\n");
	}

	fclose (fp);
}
