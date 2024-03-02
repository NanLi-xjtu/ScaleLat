#include "eam_alloy.h"

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

Setfl *setfl;

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO setfl file
------------------------------------------------------------------------- */

void allocate_alloy(char *filename)
{
  int i, j, k;
  FILE *eam_file;
  char line[1024];

  eam_file = ReadFile(filename);
  for (i = 0; i < 4; i ++) fgets(line, 1024, eam_file);
  sscanf (line, "%d", &ntypes);

  int n = ntypes;

  allocated = 1;
  AllocMem (map, n+1, int);
  for (int i = 1; i <= n; i++) map[i] = -1;
  
  AllocMem (type2frho, n+1, int);
  AllocMem2 (type2rhor, n+1, n+1, int);
  AllocMem2 (type2z2r, n+1, n+1, int);
  AllocMem2 (scale, n+1, n+1, double);

  fclose (eam_file);
}

void coeff_alloy(char *filename)
{
  int i,j;

  if (!allocated) allocate_alloy(filename);

//  if (narg != 3 + atom->ntypes)
//    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

//  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
//    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EAM setfl file
/*
  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }
*/
//  setfl = new Setfl();
  AllocMem (setfl, 1, Setfl);
//  read_file(arg[2]);
  read_file_alloy (filename);
  

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 1; i <= ntypes; i ++) map[i] = i-1;
/*
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i],setfl->elements[j]) == 0) break;
    if (j < setfl->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in EAM potential file");
  }
*/

  int n = ntypes;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      scale[i][j] = 1.0;
    }
  }

//  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

void read_file_alloy(char *filename)
{
  Setfl *file = setfl;

  // open potential file

//  int me = comm->me;
  int me = 0;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = open_potential(filename);
    if (fptr == NULL) {
      printf("Cannot open EAM potential file %s",filename);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    puts(line);
    fgets(line,MAXLINE,fptr);
    puts(line);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    n = strlen(line) + 1;
  }
//  MPI_Bcast(&n,1,MPI_INT,0,world);
//  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
//  int nwords = atom->count_words(line);
//  if (nwords != file->nelements + 1)
//    error->all(FLERR,"Incorrect element names in EAM potential file");

  char **words;
//  words = new char*[file->nelements+1];
  AllocMem (words, file->nelements+1, char*);
  int nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

//  file->elements = new char*[file->nelements];
  AllocMem (file->elements, file->nelements, char*);
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
//    file->elements[i] = new char[n];
    AllocMem (file->elements[i], n, char);
    strcpy(file->elements[i],words[i]);
  }
//  delete [] words;
  free (words);

  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    nwords = sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  if ((nwords != 5) || (file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
    printf ("Invalid EAM potential file");


  AllocMem (file->mass, file->nelements, double);

  int i, j, k;
  AllocMem2 (file->frho, file->nelements, file->nrho+1, double);
  AllocMem2 (file->rhor, file->nelements, file->nr+1, double);
  AllocMem3 (file->z2r, file->nelements, file->nelements, file->nr+1, double);

  int tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }

    if (me == 0) grab(fptr,file->nrho,&file->frho[i][1]);
    if (me == 0) grab(fptr,file->nr,&file->rhor[i][1]);
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr,file->nr,&file->z2r[i][j][1]);
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

void file2array_alloy ()
{
  int i,j,k,m,n;
  int ntypes = setfl->nelements;

  nrho = setfl->nrho;
  nr = setfl->nr;
  drho = setfl->drho;
  dr = setfl->dr;
  rhomax = (nrho-1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  AllocMem2(frho, nfrho, nrho+1, double);

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = setfl->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
//  memory->destroy(rhor);
//  memory->create(rhor,nrhor,nr+1,"pair:rhor");
  AllocMem2(rhor, nrhor, nr + 1, double);

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements+1) / 2;
//  memory->destroy(z2r);
//  memory->create(z2r,nz2r,nr+1,"pair:z2r");
  AllocMem2 (z2r, nz2r, nr + 1, double);

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}
double init_one_alloy(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setfl) {
    cutmax = 0.0;
    for (int m = 0; m < ntypes; m++)
      cutmax = MAX(cutmax,setfl[m].cut);
  } 

  cutforcesq = cutmax*cutmax;

  nmax = nMol;

  return cutmax;
}

void init_style_alloy ()
{
  // convert read-in file(s) to arrays and spline them

  file2array_alloy();
  array2spline();

  cutmax = init_one_alloy (1, 1);
}

void printfeam_alloy ()
{
  int i, j, n;
  FILE *frhoeam_file, *zream_file, *rhoream_file;

  if ((frhoeam_file = fopen ("out/KMC/eamfrho_alloy.dat", "a+")) == NULL){
    printf ("\nopen frho file error");
    getchar ();
    exit (1);
  }
  if ((zream_file = fopen ("out/KMC/eamz2r_alloy.dat", "a+")) == NULL){
    printf ("\nopen zream file error");
    getchar ();
    exit (1);
  }
  if ((rhoream_file = fopen ("out/KMC/eamrhor_alloy.dat", "a+")) == NULL){
    printf ("\nopen rhoream file error");
    getchar ();
    exit (1);
  }

  fprintf(frhoeam_file, "rho f(eV)\n");
  fprintf(rhoream_file, "r(A) rho\n");
  fprintf(zream_file, "r(A) z\n");
  for (i = 0; i < setfl->nelements; i ++) {
	fprintf (frhoeam_file, "%s\n", setfl->elements[i]);
	fprintf (rhoream_file, "%s\n", setfl->elements[i]);
	for (n = 1; n < setfl->nrho + 1;n ++){
		fprintf(frhoeam_file, "%f %f\n", setfl->drho * (n - 1), setfl->frho[i][n]);
	} 
	  
	for (n = 1; n < setfl->nr + 1;n ++){
		fprintf(rhoream_file, "%f %f\n", setfl->dr * (n - 1), setfl->rhor[i][n]);
	}
	for (j = 0; j <= i; j ++) {
		fprintf (zream_file, "%s %s\n", setfl->elements[i], setfl->elements[j]);
		for (n = 1; n < setfl->nr + 1;n ++){
			fprintf(zream_file, "%f %f\n", setfl->dr * (n - 1), setfl->z2r[i][j][n]);
		}
	}
  }

  fclose (frhoeam_file);
  fclose (zream_file);
  fclose (rhoream_file);
}
void printfspline_alloy (double ***xxxspline)
{
  double *coeff;
  double p, rsqspline, rhospline, frhop, rhop, r, z2p, z2, recip, phi, phip;
  int m, n, i, itype, jtype, num, type;
  FILE *spline_file;

  if ((spline_file = fopen ("out/KMC/eamspline.dat", "a+")) == NULL){
    printf ("\nopen spline file error");
    getchar ();
    exit (1);
  }

  num = 5000;
  AllocMem (spline, num, double);
  if (xxxspline == rhor_spline){
    itype = 1;
    for (jtype = 1; jtype <= ntypes; jtype ++) {
	    fprintf(spline_file, "r(A) rho rhop\n");
	    fprintf (spline_file,"%s\n", setfl->elements[jtype-1]);
	    for (n = 0; n < num; n++){
		rsqspline = Sqr (n * setfl[0].cut / num);
		p = sqrt(rsqspline)*rdr + 1.0;
		m = (int) (p);
		m = MIN(m,nr-1);
		p -= m;
		p = MIN(p,1.0);
		coeff = xxxspline[type2rhor[jtype][itype]][m];
		rhop = (coeff[0]*p + coeff[1])*p + coeff[2];
		spline[n] = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		fprintf (spline_file,"%e %e %e\n", sqrt(rsqspline), spline[n], rhop);
	    }
    }
  }else if (xxxspline == frho_spline){
    for (type = 1; type <= setfl->nelements; type ++) {
	    fprintf(spline_file, "rho F(eV) Fp\n");
	    fprintf (spline_file,"%s\n", setfl->elements[type-1]);
	    for (i = 0; i < num; i ++){
	      rhospline = i * setfl[0].nrho * setfl[0].drho / num; 
	      p = rhospline*rdrho + 1.0;
	      m = (int) (p);
	      m = MAX(1,MIN(m,nrho-1));
	      p -= m;
	      p = MIN(p,1.0);
	      coeff = xxxspline[type2frho[type]][m];
	      frhop = (coeff[0]*p + coeff[1])*p + coeff[2];
	      spline[i] = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	      if (rhospline > rhomax) spline[i] += frhop * (rhospline-rhomax);
	      spline[i] *= scale[1][1];
	      fprintf (spline_file,"%e %e %e\n", rhospline, spline[i], frhop);
	    }
    }
  }else if (xxxspline == z2r_spline){
	for (itype = 1; itype <= ntypes; itype ++) {
		for (jtype = 1; jtype <= itype; jtype ++) {
		    fprintf (spline_file, "r(A) z2 z2p\n");
		    fprintf (spline_file, "%s %s\n", setfl->elements[itype-1], setfl->elements[jtype-1]);
		    for (n = 0; n < num; n++){
			rsqspline = Sqr (n * setfl[0].cut / num);
			r = sqrt(rsqspline);
			p = r*rdr + 1.0;
			m = (int) (p);
			m = MIN(m,nr-1);
			p -= m;
			p = MIN(p,1.0);
			coeff = xxxspline[type2z2r[itype][jtype]][m];
			z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
			z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
			fprintf (spline_file,"%e %e %e\n", sqrt(rsqspline), z2, z2p);
		    }
		}
	}
  }

  fclose (spline_file);
}

void ComputeForcesEamPoten_alloy ()
{
  VecR drVec, DR_A;
  double fcVal;
  double *coeff;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi,fpair;
  int m, n, itype, jtype, i, j, k, j1, j2;

  AllocMem (rho, nmax, double);
  AllocMem (fp, nmax, double);

  // zero out density

  DO_MOL rho[n] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  #pragma omp parallel for private(j1, j2, i, j, drVec, DR_A, rsq, itype, jtype, k, p, m, coeff) num_threads(nthreads)
  for (n = 0; n < nebrTabLen; n ++) {
    j1 = nebrTab[2 * n];
    j2 = nebrTab[2 * n + 1];
    if (strcmp(mol[j1].elem, "vacancy") == 0 || strcmp(mol[j2].elem, "vacancy") == 0) continue;
    i = j1;
    j = j2;
    VSub (drVec, mol[j1].r, mol[j2].r);
    VWrapAll (drVec);
    VSCopy(DR_A, (lUnit * 1.e10), drVec);
    rsq = VLenSq (DR_A);
    if (rsq < cutforcesq){
      itype = jtype = 0;
      for (k = 1; k <= ntypes; k ++) {
        if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
        if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
      }
      if (itype == 0 || jtype == 0) {
        printf ("error1(eam_alloy.c):itype jtype\n");
        exit (1);
      }
      p = sqrt(rsq)*rdr + 1.0;
      m = (int) (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);
      coeff = rhor_spline[type2rhor[jtype][itype]][m];
      #pragma omp atomic
      rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      coeff = rhor_spline[type2rhor[itype][jtype]][m];
      #pragma omp atomic
      rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    }
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  frho_tot = 0.;
  #pragma omp parallel for reduction(+:frho_tot)private(itype, k, p, m, coeff, phi) num_threads(nthreads)
  for (i = 0; i < nMol; i ++){
    if (strcmp(mol[i].elem, "vacancy") == 0) continue;
    itype = 0;
    for (k = 1; k <= ntypes; k ++) {
      if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
    }
    if (itype == 0) {
        printf ("error2(eam_alloy.c):itype\n");
        exit (1);
    }
    p = rho[i]*rdrho + 1.0;
    m = (int) (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[itype]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
    phi *= scale[1][1];
    frho_tot += phi;
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  DO_MOL VZero (mol[n].ra);
  pairpoten_tot = 0.;
  uSum = 0.;
  #pragma omp parallel for reduction(+:pairpoten_tot) private(j1, j2, i, j, drVec, DR_A, rsq, itype, jtype, k, r, p, m, coeff, rhoip, rhojp, z2p, z2, recip, phi) num_threads(nthreads)
  for (n = 0; n < nebrTabLen; n ++) {
    j1 = nebrTab[2 * n];
    j2 = nebrTab[2 * n + 1];
    if (strcmp(mol[j1].elem, "vacancy") == 0 || strcmp(mol[j2].elem, "vacancy") == 0) continue;
    i = j1;
    j = j2;
    VSub (drVec, mol[j1].r, mol[j2].r);
    VWrapAll (drVec);
    VSCopy (DR_A, (lUnit * 1.e10), drVec);
    rsq = VLenSq (DR_A);
    if (rsq < cutforcesq){
        itype = jtype = 0;
        for (k = 1; k <= ntypes; k ++) {
          if (strcmp(mol[i].elem, setfl->elements[k-1]) == 0) itype = k;
          if (strcmp(mol[j].elem, setfl->elements[k-1]) == 0) jtype = k;
        }
        if (itype == 0 || jtype == 0) {
          printf ("error3(eam_alloy.c):itype jtype\n");
          exit (1);
        }
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = (int) (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip; //--eV
        pairpoten_tot += phi; //--eV
    }
  }
  uSum = frho_tot + pairpoten_tot; //--eV
//  uSum = pairpoten_tot; //--eV
//  printf("frho = %f eV, pairpoten = %f eV, uSum = %f eV\n", frho_tot, pairpoten_tot, uSum);
  uSum = uSum / (eUnit / eleChar);//dimensionless

  free (rho);
  free (fp);
} 

