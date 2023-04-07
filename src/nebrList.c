#include "nebrList.h"

int nebrListFlag = 0;
//MD parameters
double dispHi, rNebrShell, rCut;
int nMol, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;
Mol *mol;
VecI cells;
int *cellList;

void BuildNebrList ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

	rrNebr = Sqr (rCut + rNebrShell);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
	DO_MOL{
//		VSAdd (rs, mol[n].r, 0.5, region);
		VCopy (rs, mol[n].r);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;
		if (c >= (VProd(cells) + nMol) || c < 0) {
			printf ("error(nebrList.c): cellList is not big enough(surface.c),\n \
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
								VSub (dr, mol[j1].r, mol[j2].r);
								VVSub (dr, shift);
								if (VLenSq (dr) < rrNebr){
									if (nebrTabLen >= nebrTabMax) {
										printf ("error(nebrLists.c): ERR_TOO_MANY_NEBRS");
										exit (1);
									}
									nebrTab[2 * nebrTabLen] = j1;
									nebrTab[2 * nebrTabLen + 1] = j2;
									++ nebrTabLen;
								}
							}
						}
					}
				}
			}
		}
	}
}

void CheckNebrList ()
{
	int m, n, j1, j2;
	int *Nnebr, N;
	VecR dr;
	double rr, rrNebr;

	printf ("\n\n-- check neighbour list...\n\n");

	AllocMem (Nnebr, nMol, int);

	DO_MOL Nnebr[n] = 0;
	for (n = 0; n < nebrTabLen; n ++) {
		j1 = nebrTab[2 * n];
		j2 = nebrTab[2 * n + 1];
		Nnebr[j1] ++;
		Nnebr[j2] ++;
	}

	rrNebr = Sqr (rCut + rNebrShell);
	for (m = 0; m < nMol; m ++) {
		N = 0;
		for (n = 0; n < nMol; n ++) {
			if (m == n) continue;
			VSub (dr, mol[m].r, mol[n].r);
			VWrapAll (dr);
			rr = VLenSq (dr);
			if (rr < rrNebr) {
				N ++;
			}
		}
//		printf ("%d %d\n", Nnebr[m], N);
		if (Nnebr[m] != N) {
			printf ("error: nebrlist is not correct, id %d nebr %d %d\n", m, Nnebr[m], N);
			exit (1);
		}
	}

	free (Nnebr);
}

void GetNebrs (double rNebr)
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset, i, j;

	rrNebr = Sqr (rNebr);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
	DO_MOL{
//		VSAdd (rs, mol[n].r, 0.5, region);
		VCopy (rs, mol[n].r);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;
		if (c >= (VProd(cells) + nMol) || c < 0) {
			printf ("error(nebrList.c): cellList is not big enough(surface.c),\n \
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
	DO_MOL mol[n].Nnebr = 0;
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
								VSub (dr, mol[j1].r, mol[j2].r);
								VVSub (dr, shift);
								if (VLenSq (dr) < rrNebr){
									i = mol[j1].Nnebr;
									j = mol[j2].Nnebr;
									mol[j1].nebr_id[i] = j2;
									mol[j2].nebr_id[j] = j1;
									mol[j1].Nnebr ++;
									mol[j2].Nnebr ++;
									if (mol[j1].Nnebr >= nebrTabFac || mol[j2].Nnebr >= nebrTabFac) {
										printf ("error(nebrList.c): nebrTabFac %d (%d %d) is not big enough\n",\
											 nebrTabFac, mol[j1].Nnebr, mol[j2].Nnebr);
										exit (1);
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

