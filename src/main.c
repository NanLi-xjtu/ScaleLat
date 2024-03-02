#include "main.h"

char mode[128];

int main (int argc, char **argv)
{
	int time, i, j, k;
	double start, end;

	start = omp_get_wtime();
	printf("\n");
	printf("                                \\\\--//          \n");
	printf("                              (|-@ @--|)\n");
	printf("                               \\ (_)  /      \n");
	printf("                                  -         \n");
	printf("+------------o00o-------------------------------+\n");
	printf("|                                               |\n");
	printf("|                  ScaleLat                     |\n");
	printf("|        Written by Bing Xiao && Nan Li         |\n");
	printf("|                  2023-04-1                    |\n");
	printf("|  bingxiao84@xjtu.edu.cn, Bing Xiao from XJTU  |\n");
	printf("| ln906061119@stu.xjtu.edu.cn, Nan Li from XJTU |\n");
	printf("|                                               |\n");
	printf("+----------------------------------oOOo---------+\n");
	printf("\n");
	GetNameList (argc, argv);
	PrintNameList (stdout);

	Run (argc, argv);

	end = omp_get_wtime();
	time = end - start;
	printf("\ntime taken = %d h %d min %d s\n", time/3600, (time%3600)/60, time%60);
	return 0;
}

void Run (int argc, char **argv)
{
	int i;
	double T;
	char filename[128], command[128];

	srand (seed);

	switch (argc)
	{
		case 2:
			if (strcmp (argv[1], "kMC") == 0) {
				sprintf (mode, "kMC");
				system ("rm -rf ../out");
				system ("mkdir ../out");
				system ("mkdir ../out/superCell");
				system ("mkdir ../out/KMC");
				//creat super cell
				CreatSuperCell ();
				//set parameters and allocate memory
				Setup ();
				break;
			} else if (strcmp (argv[1], "restart") == 0) {
				system ("rm -rf ../out");
				system ("mkdir ../out");
				system ("mkdir ../out/KMC");
				Restart ();
			} else if (strcmp (argv[1], "analysis") == 0) {
				system ("rm -f ../out/d_profile*.dat");
				IdenPrecipitates ();
				exit (0);
			} else if (strcmp (argv[1], "MDS") == 0) {
//				error_g ();

				system ("rm -rf ../out/MDS");
				system ("mkdir ../out/MDS");
				system ("mkdir ../out/MDS/POSCAR");
				system ("cp ../xyz_to_poscar.c ../out/MDS/POSCAR");
				system ("cp ../xyz_to_poscar.sh ../out/MDS/POSCAR");
				system ("cd ../out/MDS/POSCAR && bash xyz_to_poscar.sh && cd ../../../build");

				i = 1;
				while (1) {
					sprintf (filename, "../out/MDS/POSCAR/Lattice_vectors-%d.dat", i);
					if (access(filename, 0) != 0) break;
					sprintf (command, "cp ../out/MDS/POSCAR/Lattice_vectors-%d.dat ../out/MDS/Lattice_vectors.dat", i);
					system (command);
					sprintf (command, "cp ../out/MDS/POSCAR/POSCAR-%d.dat ../out/MDS/POSCAR.dat", i);
					system (command);
					MassDensityDistribution();
					system ("cat ../out/MDS/Mass_density-1.dat >> ../out/MDS/Mass_density.movie");
					system ("cat ../out/MDS/Mass_density_A-1.dat >> ../out/MDS/Mass_density_A.movie");
					system ("cat ../out/MDS/Mass_density_B-1.dat >> ../out/MDS/Mass_density_B.movie");
					system ("cat ../out/MDS/Mass_density_AB-1.dat >> ../out/MDS/Mass_density_AB.movie");

					i ++;
				}
				exit (0);
			} else if (strcmp (argv[1], "supercell") == 0) {
				system ("rm -rf ../out");
				system ("mkdir ../out");
				system ("mkdir ../out/superCell");
				system ("mkdir ../out/KMC");
				//creat super cell
				CreatSuperCell ();
				//set parameters and allocate memory
				Setup ();
				SetupMD ();
				PrintMovie ();
				exit (0);
			} else if (strcmp (argv[1], "SRO") == 0) {
				system ("rm -f ../out/sro*.dat");
				CalSROMovie_TOTAL ();
				exit (0);
			} else if (strcmp (argv[1], "ACF") == 0) {
				system ("rm -f ../out/acf*.dat");
				CalAutocorrelationFunction ();
				exit (0);
			} else if (strcmp (argv[1], "CE") == 0) {
				system ("rm -rf ../out && mkdir ../out");
				system ("mkdir ../out/CE_benchmark");
				ReadMat ();
				ObtainSym (); //obtain symmetry points group
				GoalClusterExtract ();
				exit (0);
			} else if (strcmp (argv[1], "MAP") == 0) {
				sprintf (mode, "MAP");
				system ("rm -rf ../out && mkdir ../out");
				system ("mkdir ../out/CE_benchmark");
				system ("mkdir ../out/CE_trial");
				system ("mkdir ../out/superCell");
				ReadMat ();
				ObtainSym (); //obtain symmetry points group
				CreatSuperCell ();//creat super cell
				system ("cp ../out/superCell/Supercell.xyz ../in/trial.xyz");
				system ("rm -rf ../out/superCell");
				Mapping ();
				exit (0);
			} else {
				printf ("error(mc.c): input format\n");
				exit (1);
			}
			break;
		default:
			printf ("error(mc.c): input format\n");
			exit (1);

	}

	SetupMD ();

	/***********************Monte Carlo calculation***********************/
	printf ("\n-- Monte Carlo calculation\n");
	//initialize
	T = T_max; //--K
	printf ("temperature: %.4f K\n", T);

	//main loop
	if (strcmp (MC_method, "vacancy") == 0) MonteCarloDoubleT (T); //--K
	else if (strcmp (MC_method, "atom") == 0) MonteCarloExchangeAtoms (T);
	else if (strcmp (MC_method, "mix") == 0) MonteCarloMixed (T);
	else {
		printf ("error(mc.c): MC_method format\n");
		exit  (1);
	}

}
