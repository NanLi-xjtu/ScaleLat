// The code is written by Bing Xiao, my current address is: 
// Depertment of Earth Sciences
// University College London
// Gower Street, London, W1CE 6BT, United Kingdom
// b.xiao@ucl.ac.uk
/* This program calculates the mass density distribution for 3-D crystal structures or
any atomic structures with periodic boundary conditions */
/*
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <iostream>
*/

#include "Mass_density_distribution.h"

#define MAX_MDS 10000000
#define NMAX 10000000  // Atomic coordinates in supercell
#define delta 2.4
#define PI 3.1415926508
#define Dalton 1.660539E-24
#define Vcm 1.0E-24 
int MassDensityDistribution()
{
      int i,j,k,l,t,r;
      int I,J,K,L;
      int R_1, R_2, R_3;
      int tmp_1, tmp_2, tmp_3;
      int frame; // Total number of frames for the calculations
      int methods;
      int options;
      char finame[64], finame1[64], finame2[64], finame3[64], finame4[64]; // File names
      char finame5[64], finame6[64], finame7[64],finame8[64]; // File names
      char files[128], files1[128];
      printf("/*----------------------------------------------------------------*/\n");
      printf("/*      THIS PROGRAM IS WRITTEN BY BING XIAO                      */\n");
      printf("/*-----------------CURRENT ADDRESS--------------------------------*/\n");
      printf("/*           Department of Earth Sciences                         */\n");
      printf("/*            University College London                           */\n");
      printf("/*      Gower Street, London WC1E 6BT, United Kingdom             */\n");
      printf("/*                  b.xiao@ucl.ac.uk                              */\n");
      printf("/*     READING XDATCAR FILES AND CALCULATING Z-Porfile            */\n");
      printf("/*         Convert XDATCAR to xyz.dat fomrmat                     */\n");
      printf("/*            Center of Mass for each snapshot                    */\n");
      printf("/*             3-D Mass Density Distribution                      */\n");
      printf("/*           Isovalue 3-D Contour Distribution                    */\n"); 
      printf("/*----------------------------------------------------------------*/\n");
      printf("/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/\n");
      printf("/*@  Combining Multiple XDATCARs to a Single XDATCAR File        @*/\n");
      printf("/*@  The First Seven HeadLines MUST BE REMOVED except the first  @*/\n");
      printf("/*@  XDATCAR file,otherwise,the program does not run properly.   @*/\n");
      printf("/*@  Using Linux Command: sed -i '1,7d' XDATCAR_xxx              @*/\n");
      printf("/*@  Then combine all XDATCARs: cat XDATCAR_xxx* > XDATCAR       @*/\n");
      printf("/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/\n"); 
      printf("IN THE CURRENT IMPLEMENTATION, CALCULATIONS ARE LIMITED TO SIX\n");
      printf("ELEMENTS.\n");
// Define the lattic vectors
      double a_1, a_2, a_3; //the lattice vector in the first direction [100] and its conponents on three principal axies.
      double b_1, b_2, b_3; // [010]
      double c_1, c_2, c_3; // [001]
      float *Latx, *Laty, *Latz; // The lattice vectors stored in the input file
      double a, b, c; // Cell constants
      double start_c, final_c;
      double Grid_x, Grid_y, Grid_z; // The discretized grid in all directions
      double step_x, step_y, step_z; // The step of discretized lattice volume
      double V; // Cell volume
      double V_e; // Volume element
// Define the atomic number and atomic postions in the crystal structure
      int *Name; // Atomic number
      int atomA, atomB, atomC; // atomic number
      int atomD, atomE, atomF;
      float *fracx, *fracy, *fracz;

      AllocMem (Latx, MAX_MDS, float);
      AllocMem (Laty, MAX_MDS, float);
      AllocMem (Latz, MAX_MDS, float);
      AllocMem (Name, MAX_MDS, int);
      AllocMem (fracx, MAX_MDS, float);
      AllocMem (fracy, MAX_MDS, float);
      AllocMem (fracz, MAX_MDS, float);

      double Atonum_1, Atonum_2, Atonum_3; // Atomic weight of each species
      double Atonum_4, Atonum_5, Atonum_6;
      double mass_A, mass_B, mass_C, mass_t; // Mass of each element in the cell and total mass
      double mass_D, mass_E, mass_F;
      double massA_cal, massB_cal, massC_cal; // calculated by smearing
      double massD_cal, massE_cal, massF_cal;
      double mass_cal, mass_calt; // Calculated mass in the cell after smearing
      double normal; // The normalizing constant 
      printf("Atomic numbers for Si, Mg, O are 14, 12 and 8\n");
      printf("Atomic masses are 28.086, 24.305 and 15.999\n");
      printf("Input the parameters for all six elements below.\n");
      printf("------------------------------------------------------------------\n");
      printf("If the elements are less than 6, please use 0 for their parameters.\n");
      printf("Atomic weight for first three elements:(15.999 28.086 24.305)\n");
//      scanf("%lf %lf %lf",&Atonum_1, &Atonum_2, &Atonum_3);
      Atonum_1 = 55.84;
      Atonum_2 = 52.00;
      Atonum_3 = 0.;
      printf("Atomic weight for last three elements:(1.000 19.000 56.000)\n");
//      scanf("%lf %lf %lf",&Atonum_4, &Atonum_5, &Atonum_6);
      Atonum_4 = 0.;
      Atonum_5 = 0.;
      Atonum_6 = 0.;
      printf("------------------------------------------------------------------\n");
      printf("Enter atomic number for first three elements:(8 14 0)\n");
//      scanf("%d %d %d",&atomA, &atomB, &atomC);
      atomA = 26;
      atomB = 24;
      atomC = 0;
      printf("Enter atomic number for last three elements:(1 19 26)\n");
//      scanf("%d %d %d",&atomD, &atomE, &atomF);
      atomD = 0;
      atomE = 0;
      atomF = 0;
      printf("------------------------------------------------------------------\n");
// Define the supercell size
      int *D_1;
      int C_x, C_y, C_z; // Supercell dimension
      int Cell_x, Cell_y, Cell_z;
      float *Coordx, *Coordy, *Coordz;

      AllocMem (D_1, NMAX, int);
      AllocMem (Coordx, NMAX, float);
      AllocMem (Coordy, NMAX, float);
      AllocMem (Coordz, NMAX, float);

      float X_1, Y_1, Z_1;
      float G_1x, G_1y, G_1z;
// Define the cubic grids in the unit cell
      float X_cell, Y_cell, Z_cell;
      float r_dis;
// Auxiliary vribales and functions
      double x, y, z;
      double x_1, y_1, z_1;
      float Zcoord; // Z ccordinate
      float f(double, double, double);
      float g(double, double); // Mass smearing function
      double *rou_A, *rou_B, *rou_C, *rou; // Mass density at r for each species and sum of them
      double *rou_D, *rou_E, *rou_F;

      AllocMem (rou, NMAX, double);
      AllocMem (rou_A, NMAX, double);
      AllocMem (rou_B, NMAX, double);
      AllocMem (rou_C, NMAX, double);
      AllocMem (rou_D, NMAX, double);
      AllocMem (rou_E, NMAX, double);
      AllocMem (rou_F, NMAX, double);

      double isovalue1, isovalue2, isovalue3;
      double range_1, range_2, range_3;
      double range_4; // range_5;
      FILE *input;
      FILE *input1;
      FILE *output1;      
      FILE *output2;
      FILE *output3;
      FILE *output4;
      FILE *output5;
      FILE *output7;
      FILE *output8;
      input=fopen("../out/MDS/Lattice_vectors.dat","r");
     /* READ LATTICE VECTORS */
      i=0;
      while(fscanf(input,"%f %f %f",&Latx[i],&Laty[i],&Latz[i])!=EOF){
                              i=i+1;
                              }
      a_1=Latx[0]; b_1=Latx[1]; c_1=Latx[2];
      a_2=Laty[0]; b_2=Laty[1]; c_2=Laty[2];
      a_3=Latz[0]; b_3=Latz[1]; c_3=Latz[2];
      a=f(a_1,a_2,a_3);
      b=f(b_1,b_2,b_3);
      c=f(c_1,c_2,c_3);
      V=a*b*c;
      fclose(input);
      /* CHECK THE INPUTS */
      printf("The lattice vectors are\n");
      printf("%lf\t%lf\t%lf\n",a_1,a_2,a_3);
      printf("%lf\t%lf\t%lf\n",b_1,b_2,b_3);
      printf("%lf\t%lf\t%lf\n",c_1,c_2,c_3);
      printf("The lattice constants are\n");
      printf("a (Angstrom)=%lf\n",a);
      printf("b (Angstrom)=%lf\n",b);
      printf("c (Angstrom)=%lf\n",c);
      printf("Cell Volume is=%f\n",V);
      printf("Atomic weigh of species A is=%lf\tB is=%lf\t C is=%lf\n",Atonum_1,Atonum_2, Atonum_3);
      printf("Atomic weigh of species D is=%lf\tE is=%lf\t F is=%lf\n",Atonum_4,Atonum_5, Atonum_6);
     /* READ ATOMIC POSITIONS IN SINGLE CELL */
      printf("****************************************************************\n");
      printf("THIS IS A TESTING PROGRAM FOR SINGLE MD FRAME!\n");
      printf("THE TOTAL FRAME NUMBER SHOULD BY ONE\n");
      printf("FOR MULTIPLE MD FRAMES, USING ANOTHER PROGRAM\n");
      printf("****************************************************************\n");
      printf("Please indicate the total number of frame:(MUST BE 1)\n");
//      scanf("%d",&frame);
      frame = 1;
      printf("Enter grid points in three directions (10.0 20.0 100.0)\n");
      printf("Decease this step size will significanly increase the total number of data points computed\n");
      printf("Thus increasing the computational time\n");
//      scanf("%lf %lf %lf", &Grid_x, &Grid_y, &Grid_z);
      Grid_x = Grid_y = Grid_z = 100;
      printf("Three isovalue surfaces are calculated, user need specify three isovalues below.\n");
      printf("If you wish to calcfulate only one isovalue surface, use the same numbers for other two.\n");
      printf("Giving a value for the first isovalue surface(Liquid):\n");
//      scanf("%lf",&isovalue1);
      isovalue1 = 1.;
      printf("The isovalue width:\n");
//      scanf("%lf",&range_1);
      range_1 = 1.;
      printf("Fluctuation amplitude:\n");
//      scanf("%lf",&range_4);
      range_4 = 1.;
      printf("Second isovalue(Vapor):\n");
//      scanf("%lf",&isovalue2);
      isovalue2 = 1.;
      printf("The isovalue width:\n");
//      scanf("%lf",&range_2);
      range_2 = 1.;
      printf("Third isovalue:\n");
//      scanf("%lf",&isovalue3);
      isovalue3 = 1.;
      printf("The isovalue width:\n");
//      scanf("%lf",&range_3);
      range_3 = 1.;
      printf("------------------------------------------------------------------\n");
      printf("c lattice constant is assumed to be the longest one.\n");
      printf("You may specify the c range in the calculation.\n");
      printf("Two options: Default range(1); User defined range (2)\n");
      printf("Enter you choice:1 or 2\n");
//      scanf("%d",&options);
      options = 1;
      step_x=a/Grid_x;
      step_y=b/Grid_y;
      if(options==1){
                     start_c=0.00;
                     final_c=c;
                     step_z=c/Grid_z;
                     }
      else if(options==2){
                          printf("You must specify the range for c lattice constant.\n");
                          printf("Enter c range below:(0.0000 %lf)\n",c);
                          scanf("%lf %lf",&start_c,&final_c);
                          step_z=fabs((start_c-final_c)/Grid_z);
                          }
      printf("The unit is Angstrom\n");
      printf("Step size in x direction is=%lf\n",step_x);
      printf("Step size in y direction is=%lf\n",step_y);
      printf("Step size in z direction is=%lf\n",step_z);
      printf("The grid in X direction is=%lf\nThe grid in Y direction is=%lf\n",Grid_x+1.0, Grid_y+1.0);
      printf("The grid in Z direction is=%lf\n",Grid_z+1.0);
      printf("RUNNING, PLEASE WAIT\n");
/*-------------------------GENERATING THE CUBIC GRID FOR UNIT CELL--------------*/
      FILE *output6; // Save the cubic grid coordinates
      output6=fopen("../out/MDS/Cubic_grid.dat","w");
      for(X_cell=0.00;X_cell<=a;X_cell=X_cell+step_x){
                                                      for(Y_cell=0.00;Y_cell<=b;Y_cell=Y_cell+step_y){
                                                                                                      for(Z_cell=start_c;Z_cell<=final_c;Z_cell=Z_cell+step_z){
                                                                                                                                                      fprintf(output6,"%f\t%f\t%f\n",X_cell, Y_cell, Z_cell);
                                                                                                                                                      }
                                                                                                             }
                                                           }
      fclose(output6);
// Check the total number of points in the cubic_gird.dat file
      float spx, spy, spz; // dummy varibales
      int spnum, rowsp;
      float *super_x, *super_y, *super_z;
      output6=fopen("../out/MDS/Cubic_grid.dat","r");
/* Find how many rows in the xyz file */
      spnum=0;
      while(!feof(output6)){
                            fscanf(output6,"%f %f %f",&spx, &spy, &spz);
                            spnum++;
                            }
      rowsp=spnum-1;
      fclose(output6);
      printf("Total number of gird points in the unit cell is %d\n",rowsp);
// Allocate the memoy to the pointer
      super_x=(float*)malloc(rowsp*sizeof(float));
      super_y=(float*)malloc(rowsp*sizeof(float));
      super_z=(float*)malloc(rowsp*sizeof(float));
      if(super_x==NULL||super_y==NULL||super_z==NULL){
                                                   printf("Error-could not allocate an array.\n");
                                                   }
      printf("Starting to read to cubic grid.\n");
      output6=fopen("../out/MDS/Cubic_grid.dat","r");
      for(i=0;i<rowsp;i++){
                          fscanf(output6,"%f %f %f",&(super_x[i]),&(super_y[i]),&(super_z[i]));
                          }
      fclose(output6);
      printf("All coordinates in the cubic grids are stored in the array.\n");
//------------------------------------------------------------------------------
/*---------------------------MASS DENSITY DISTRIBUTION--------------------------*/
      for(L=1;L<=frame;L++){
                         sprintf(finame,"../out/MDS/POSCAR.dat");
                         input1=fopen(finame,"r");
                         if(input1==NULL){
                                          printf("Can't open %s\n",finame);
                                          exit(1);
                                          }
                         j=0;
                         while(fscanf(input1,"%d %f %f %f",&Name[j],&fracx[j],&fracy[j],&fracz[j])!=EOF){
                                j=j+1;
                                }
      R_1=j;
      fclose(input1);
      // Check the coordinates again
      printf("The atomic postions are given blow for the current cell\n");
      mass_A=0.0;
      mass_B=0.0;
      mass_C=0.0; mass_D=0.0; mass_E=0.0; mass_F=0.0;
      for(i=0;i<R_1;i++){
                         if(Name[i]==atomA){
                                            mass_A+=Atonum_1;
                                            }
                         else if(Name[i]==atomB){
                                                 mass_B+=Atonum_2;
                                                 }
                         else if(Name[i]==atomC){
                                                 mass_C+=Atonum_3;
                                                 }
                         else if(Name[i]==atomD){
                                                 mass_D+=Atonum_4;
                                                 }
                         else if(Name[i]==atomE){
                                                 mass_E+=Atonum_5;
                                                 }
                         else if(Name[i]==atomF){
                                                 mass_F+=Atonum_6;
                                                 }
//                         printf("%d\t%f\t%f\t%f\n",Name[i],fracx[i],fracy[i],fracz[i]);
                         }
      mass_t=mass_A+mass_B+mass_C+mass_D+mass_E+mass_F;
      printf("The total number of atoms is=%d\n",R_1);
      printf("Total mass of the cell is(g/mol)=%lf\n",mass_t);
/* Define the supercell dimensions and computing the cooredinates for all atoms in the cell */
      printf("Two methods are implemented in this program to calculate mass density.\n");
      printf("Method 1: Supercell method (3x3x3 supercell will be generated to treat PBC).\n");
      printf("Method 2: Unit cell method (1x1x1 cell is used with an advanced treatment for PBC).\n");
      printf("Notice: Method 2 is 27 times faster than method 1.\n");
      printf("Choose your method: 1 or 2\n");
//      scanf("%d",&methods);
      methods = 2.;
      FILE *output;
      sprintf(finame1,"../out/MDS/Supercell-%d.dat",L);
      output=fopen(finame1,"w");
      if(methods==1){
                     printf("Input the dimensional of the supercell\n");
                     printf("Enter the three integers as 1 2 3\n");
                     printf("For the purpose of this program, please input 1 1 1\n");
                     scanf("%d %d %d",&Cell_x, &Cell_y, &Cell_z);
      for(i=0;i<R_1;i++){
                         X_1=sqrt(pow(fracx[i]*a_1,2.0)+pow(fracy[i]*a_2,2.0)+pow(fracz[i]*a_3,2.0));
                         Y_1=sqrt(pow(fracx[i]*b_1,2.0)+pow(fracy[i]*b_2,2.0)+pow(fracz[i]*b_3,2.0));
                         Z_1=sqrt(pow(fracx[i]*c_1,2.0)+pow(fracy[i]*c_2,2.0)+pow(fracz[i]*c_3,2.0));
                         R_2=Name[i];
                         for(j=-Cell_x;j<Cell_x+1;j++){
                                                      for(k=-Cell_y;k<Cell_y+1;k++){
                                                                                    for(l=-Cell_z;l<Cell_z+1;l++){
                                                                                                                  G_1x=X_1+j*sqrt(pow(a_1,2.0)+pow(a_2,2.0)+pow(a_3,2.0))+0.0*k+0.0*l;
                                                                                                                  G_1y=Y_1+k*sqrt(pow(b_1,2.0)+pow(b_2,2.0)+pow(b_3,2.0))+0.0*j+0.0*l;
                                                                                                                  G_1z=Z_1+l*sqrt(pow(c_1,2.0)+pow(c_2,2.0)+pow(c_3,2.0))+0.0*k+0.0*j;
                                                                                                               //   printf("The coordinates of the first species in supercell are\n");
                                                                                                                  printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\n", R_2,j,k,l,G_1x,G_1y,G_1z);
                                                                                                                  fprintf(output,"%d\t%f\t%f\t%f\n",R_2,G_1x,G_1y,G_1z);
                                                                                                                  }
                                                                                    }
                                                      }
                         }
                     }
      else if(methods==2){
                          Cell_x=0; Cell_y=0; Cell_z=0;
                          printf("Only the unit cell is used in the calculations.\n");
                          for(i=0;i<R_1;i++){
                                             X_1=sqrt(pow(fracx[i]*a_1,2.0)+pow(fracy[i]*a_2,2.0)+pow(fracz[i]*a_3,2.0));
                                             Y_1=sqrt(pow(fracx[i]*b_1,2.0)+pow(fracy[i]*b_2,2.0)+pow(fracz[i]*b_3,2.0));
                                             Z_1=sqrt(pow(fracx[i]*c_1,2.0)+pow(fracy[i]*c_2,2.0)+pow(fracz[i]*c_3,2.0));
                                             R_2=Name[i];
//                                             printf("%d\t%f\t%f\t%f\n", R_2,X_1,Y_1,Z_1);
                                             fprintf(output,"%d\t%f\t%f\t%f\n",R_2,X_1,Y_1,Z_1);
                                             }
                          }
      fclose(output);
      printf("The coordinates are stored in file\n");
/* ACCOUNT TOTAL NUMBER OF ATOMS IN THE SUPERCELL */
      sprintf(finame2,"../out/MDS/Supercell-%d.dat",L);
      output=fopen(finame2,"r");
      J=0;
      while(fscanf(output,"%d %f %f %f",&D_1[J],&Coordx[J],&Coordy[J],&Coordz[J])!=EOF){
                                J=J+1;
                                }
      R_3=J;
      fclose(output);
      printf("Number of atoms in the single cell=%d\n",R_1);
      printf("The total number of atoms in the supercell is=%d\n",R_3);
//      printf("The dimension of the supercell is 3 X 3 X 3\n");
// COMPUTING THE MASS DENSITY FIELD 
      sprintf(finame3,"../out/MDS/Mass-%d.dat",L);
      sprintf(finame8,"../out/MDS/Mass_species-%d.dat",L);
      sprintf(files,"../out/MDS/Mass-%d-tmp.dat",L);
      sprintf(files1,"../out/MDS/Mass_species-%d-tmp.dat",L);
      output1=fopen(finame3,"w");
      output5=fopen(finame8,"w");
      output7=fopen(files,"w");
      output8=fopen(files1,"w");
      mass_calt=0.00;
      fprintf(output1,"X\tY\tZ\tDensity(g/cm3)\n");
      fprintf(output5,"X\tY\tZ\trou_A\trou_B\trou_C\trou_D\trou_E\trou_F\n");
//      for(X_cell=0.00;X_cell<a+step_x;X_cell=X_cell+step_x){
                                                    //  for(Y_cell=0.00;Y_cell<b+step_y;Y_cell=Y_cell+step_y){
                                                                                                         //    for(Z_cell=0.00;Z_cell<c+step_z;Z_cell=Z_cell+step_z){
/*------------------------------------------------------------------------------*/
     int step = 0;
     if(methods==1){
//                  #pragma omp parallel for reduction(+:mass_calt)private(massA_cal, massB_cal, massC_cal, massD_cal, massE_cal, massF_cal, K, x_1, y_1, z_1, r_dis, V_e, mass_cal) num_threads(nthreads)
                    for(i=0;i<rowsp;i++){
                                         rou_A[i]=0.00; rou_D[i]=0.00;
                                         rou_B[i]=0.00; rou_E[i]=0.00;
                                         rou_C[i]=0.00; rou_F[i]=0.00;
                                         massA_cal=0.00; massD_cal=0.00;
                                         massB_cal=0.00; massE_cal=0.00;
                                         massC_cal=0.00; massF_cal=0.00;
                                         for(K=0;K<R_3;K++){
                                                            x_1=super_x[i]-Coordx[K];
                                                            y_1=super_y[i]-Coordy[K];
                                                            z_1=super_z[i]-Coordz[K];
                                                            r_dis=f(x_1,y_1,z_1);
//                                                            if (r_dis > 10) break;
                                                            V_e=step_x*step_y*step_z;
                                                            if(D_1[K]==atomA){
                                                                              rou_A[i]+=(Atonum_1*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                              massA_cal+=Atonum_1*g(delta,r_dis);
                                                                              }
                                                            else if(D_1[K]==atomB){
                                                                                   rou_B[i]+=(Atonum_2*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massB_cal+=Atonum_2*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomC){
                                                                                   rou_C[i]+=(Atonum_3*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massC_cal+=Atonum_3*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomD){
                                                                                   rou_D[i]+=(Atonum_4*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massD_cal+=Atonum_4*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomE){
                                                                                   rou_E[i]+=(Atonum_5*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massE_cal+=Atonum_5*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomF){
                                                                                   rou_F[i]+=(Atonum_6*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massF_cal+=Atonum_6*g(delta,r_dis);
                                                                                   }
                                                            }
                                                            mass_cal=massA_cal+massB_cal+massC_cal+massD_cal+massE_cal+massF_cal;
                                                            rou[i]=(rou_A[i]+rou_B[i]+rou_C[i]+rou_D[i]+rou_E[i]+rou_F[i]); //0.0398124;
                                                            mass_calt+=mass_cal;
                                                            printf("%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);
                                                            
                                         }
                    for(i=0;i<rowsp;i++){
                                                            fprintf(output5,"%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",super_x[i],super_y[i],super_z[i],rou_A[i],rou_B[i],rou_C[i],rou_D[i],rou_E[i],rou_F[i]);
                                                            fprintf(output1,"%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);
                                                            fprintf(output8,"%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",super_x[i],super_y[i],super_z[i],rou_A[i],rou_B[i],rou_C[i],rou_D[i],rou_E[i],rou_F[i]);
                                                            fprintf(output7,"%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);
                                        }
                    }
     else if(methods==2){
                       #pragma omp parallel for reduction(+:mass_calt)private(massA_cal, massB_cal, massC_cal, massD_cal, massE_cal, massF_cal, K, x_1, y_1, z_1, r_dis, V_e, mass_cal) num_threads(nthreads)
                         for(i=0;i<rowsp;i++){
                                         rou_A[i]=0.00; rou_D[i]=0.00;
                                         rou_B[i]=0.00; rou_E[i]=0.00;
                                         rou_C[i]=0.00; rou_F[i]=0.00;
                                         massA_cal=0.00; massD_cal=0.00;
                                         massB_cal=0.00; massE_cal=0.00;
                                         massC_cal=0.00; massF_cal=0.00;
                                        // fprintf(output1,"X\tY\tZ\tDensity(g/cm3)\n");
                                        // fprintf(output5,"X\tY\tZ\trou_A\trou_B\trou_C\trou_D\trou_E\trou_F\n");
                                         for(K=0;K<R_3;K++){
                                                            x_1=super_x[i]-Coordx[K];
                                                            y_1=super_y[i]-Coordy[K];
                                                            z_1=super_z[i]-Coordz[K];
                                                            if(x_1 < -a/2.0){
                                                                             x_1=x_1+a;
                                                                             }
                                                            else if(x_1 > a/2.0){
                                                                                 x_1=x_1-a;
                                                                                 }
                                                            if(y_1 < -b/2.0){
                                                                             y_1=y_1+b;
                                                                            }
                                                            else if(y_1 > b/2.0){
                                                                                 y_1=y_1-b;
                                                                                }
                                                            if(z_1 < -c/2.0){
                                                                             z_1=z_1+c;
                                                                             }
                                                            else if(z_1 > c/2.0){
                                                                                 z_1=z_1-c;
                                                                                 }  
                                                            r_dis=f(x_1,y_1,z_1);
                                                            if (r_dis > 10.) continue;
                                                            V_e=step_x*step_y*step_z;
                                                            if(D_1[K]==atomA){
                                                                              rou_A[i]+=(Atonum_1*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                              massA_cal+=Atonum_1*g(delta,r_dis);
                                                                              }
                                                            else if(D_1[K]==atomB){
                                                                                   rou_B[i]+=(Atonum_2*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massB_cal+=Atonum_2*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomC){
                                                                                   rou_C[i]+=(Atonum_3*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massC_cal+=Atonum_3*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomD){
                                                                                   rou_D[i]+=(Atonum_4*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massD_cal+=Atonum_4*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomE){
                                                                                   rou_E[i]+=(Atonum_5*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massE_cal+=Atonum_5*g(delta,r_dis);
                                                                                   }
                                                            else if(D_1[K]==atomF){
                                                                                   rou_F[i]+=(Atonum_6*Dalton/(1.0*Vcm))*g(delta,r_dis);
                                                                                   massF_cal+=Atonum_6*g(delta,r_dis);
                                                                                   }
                                                            }
                                                            mass_cal=massA_cal+massB_cal+massC_cal+massD_cal+massE_cal+massF_cal;
                                                            rou[i]=(rou_A[i]+rou_B[i]+rou_C[i]+rou_D[i]+rou_E[i]+rou_F[i]); //0.0398124;
                                                            mass_calt+=mass_cal;
//                                                            printf("%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);
                                                            if (i % (int)(rowsp/1000.) == 0) printf("%.1f%%\t", step++ / 10.);

                                         }

			for(i=0;i<rowsp;i++) {
 				fprintf(output5,"%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",super_x[i],super_y[i],super_z[i],rou_A[i],rou_B[i],rou_C[i],rou_D[i],rou_E[i],rou_F[i]);
				fprintf(output1,"%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);
				fprintf(output8,"%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",super_x[i],super_y[i],super_z[i],rou_A[i],rou_B[i],rou_C[i],rou_D[i],rou_E[i],rou_F[i]);
				fprintf(output7,"%f\t%f\t%f\t%lf\n",super_x[i],super_y[i],super_z[i],rou[i]);

			}

                    }


      normal=mass_t/mass_calt;
      fclose(output1);
      fclose(output5);
      fclose(output7); fclose(output8);
      printf("The normalization factor is=%f\n",normal);
// CALCULATE THE DENSITY AS A FUNCTION OF Z
      float rouz,rouz_A,rouz_B,rouz_C;
      float rouz_D, rouz_E, rouz_F;
      float weight;
      float dummy, dummy1, dummy2, dummy3;
      int rown; // Total number of rows in the data
      float *Xdata,*Ydata,*Zdata,*Wdata;
      float *XXdata, *YYdata, *ZZdata, *WAdata, *WBdata, *WCdata;
      float *WDdata, *WEdata, *WFdata;
      float gds; // Gibbs dividing surface 
      char file1[128], file2[128];
      char file3[128], file4[128];
      FILE *output9;
      FILE *output10;
      FILE *output11;
      FILE *output12;
      sprintf(finame4,"../out/MDS/density_isovalue1-%d.dat",L);
      sprintf(finame5,"../out/MDS/Mass-%d-tmp.dat",L);
      sprintf(finame6,"../out/MDS/Z_density-%d.dat",L);
      sprintf(finame7,"../out/MDS/Mass_density-%d.dat",L);
      sprintf(file1,"../out/MDS/density_isovalue2-%d.dat",L);
      sprintf(file2,"../out/MDS/density_isovalue3-%d.dat",L);
      sprintf(file3,"../out/MDS/density_liquid.dat");
      sprintf(file4,"../out/MDS/density_vapor.dat");
      output7=fopen(finame5,"r");
      output2=fopen(finame4,"w");
      output3=fopen(finame6,"w");
      output4=fopen(finame7,"w");
      output9=fopen(file1,"w");
      output10=fopen(file2,"w");
      output11=fopen(file3,"w");
      output12=fopen(file4,"w");
      // Determine the total number of the rows in the input file
      K=0;
      while(!feof(output7)){
                            fscanf(output7,"%f %f %f %f",&dummy, &dummy1, &dummy2, &dummy3);
                            K++;
                            }
      rown=K-1; 
      fclose(output7);
      // Allocate the dynamical memory for the arrays
      Xdata=(float*)malloc(rown*sizeof(float));
      Ydata=(float*)malloc(rown*sizeof(float));
      Zdata=(float*)malloc(rown*sizeof(float));
      Wdata=(float*)malloc(rown*sizeof(float));
      XXdata=(float*)malloc(rown*sizeof(float));
      YYdata=(float*)malloc(rown*sizeof(float));
      ZZdata=(float*)malloc(rown*sizeof(float));
      WAdata=(float*)malloc(rown*sizeof(float));
      WBdata=(float*)malloc(rown*sizeof(float));
      WCdata=(float*)malloc(rown*sizeof(float));
      WDdata=(float*)malloc(rown*sizeof(float));
      WEdata=(float*)malloc(rown*sizeof(float));
      WFdata=(float*)malloc(rown*sizeof(float));
      if(Xdata==NULL&&Ydata==NULL&&Zdata==NULL&&Wdata==NULL){
                                                             printf("Error-could not allocate an array.\n");
                                                             }
      output7=fopen(finame5,"r");
      for(I=0;I<rown;I++){
                          fscanf(output7,"%f %f %f %f",&(Xdata[I]),&(Ydata[I]),&(Zdata[I]),&(Wdata[I]));
                          }
      if(XXdata==NULL||YYdata==NULL||ZZdata==NULL||WAdata==NULL||WBdata==NULL||WCdata==NULL){
                                                                                              printf("Error-could not allocate an array.\n");
                                                                                              }
      output8=fopen(files1,"r");
      for(I=0;I<rown;I++){
                          fscanf(output8,"%f %f %f %f %f %f %f %f %f",&(XXdata[I]),&(YYdata[I]),&(ZZdata[I]),&(WAdata[I]),&(WBdata[I]),&(WCdata[I]),&(WDdata[I]),&(WEdata[I]),&(WFdata[I]));
                          }
      printf("Total number of row=%d\n",K);
      printf("3-D grid is=%f\t%f\t%f\n",Grid_x+1.0,Grid_y+1.0,Grid_z+1.0);
      gds = (isovalue1+isovalue2)/2.0;
      fprintf (output4, "%d\nSolutionReader properties=pos:R:3:rho:R:1\n", rown);

	FILE *OUTPUT;
	char filename[128];
	int Na, Nb, Nab;
	Na = Nb = Nab = 0;
      for(t=0;t<rown;t++){
                          fprintf(output4,"%f\t%f\t%f\t%f\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
				if (Wdata[t] < 7.408) Na ++;
				else if (Wdata[t] > 7.732) Nb ++;
				else Nab ++;
                          if(Wdata[t]>=isovalue1-range_1&&Wdata[t]<=isovalue1+range_1){
                                                                      fprintf(output2,"%f\t%f\t%f\t%lf\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
                                                                      }
                          if(Wdata[t]>=isovalue2-range_2&&Wdata[t]<=isovalue2+range_2){
                                                                      fprintf(output9,"%f\t%f\t%f\t%lf\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
                                                                      }
                          if(Wdata[t]>=isovalue3-range_3&&Wdata[t]<=isovalue3+range_3){
                                                                      fprintf(output10,"%f\t%f\t%f\t%lf\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
                                                                      }
                          if(Wdata[t]>=gds){
                                            fprintf(output11,"%f\t%f\t%f\t%lf\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
                                            }
                          if(Wdata[t]<=gds){
                                            fprintf(output12,"%f\t%f\t%f\t%lf\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
                                            }
                          }



			
			sprintf (filename, "../out/MDS/Mass_density_A-1.dat");
			OUTPUT = fopen (filename, "w");
			if (Na > 0) fprintf (OUTPUT, "%d\nSolutionReader properties=pos:R:3:rho:R:1\n", Na);
			for(t=0;t<rown;t++) {	
				if (Wdata[t] < 7.408) fprintf(OUTPUT,"%f\t%f\t%f\t%f\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
			}
			fclose (OUTPUT);
			sprintf (filename, "../out/MDS/Mass_density_B-1.dat");
			OUTPUT = fopen (filename, "w");
			if (Nb > 0) fprintf (OUTPUT, "%d\nSolutionReader properties=pos:R:3:rho:R:1\n", Nb);
			for(t=0;t<rown;t++) {	
				if (Wdata[t] > 7.732) fprintf(OUTPUT,"%f\t%f\t%f\t%f\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
			}
			fclose (OUTPUT);
			sprintf (filename, "../out/MDS/Mass_density_AB-1.dat");
			OUTPUT = fopen (filename, "w");
			if (Nab > 0) fprintf (OUTPUT, "%d\nSolutionReader properties=pos:R:3:rho:R:1\n", Nab);
			for(t=0;t<rown;t++) {	
				if (Wdata[t] >= 7.408 && Wdata[t] <= 7.732) fprintf(OUTPUT,"%f\t%f\t%f\t%f\n",Xdata[t],Ydata[t],Zdata[t],Wdata[t]);
			}
			fclose (OUTPUT);

			sprintf (filename, "../out/MDS/volume.dat");
			OUTPUT = fopen (filename, "a+");
			fprintf (OUTPUT, "%f %f\n", (double)Na/(double)rown, (double)Nb/(double)rown);
			fclose (OUTPUT);



      for(Zcoord=start_c;Zcoord<final_c+step_z;Zcoord=Zcoord+step_z){
                                                           rouz=0.0;
                                                           rouz_A=0.0; rouz_D=0.0;
                                                           rouz_B=0.0; rouz_E=0.0;
                                                           rouz_C=0.0; rouz_F=0.0;
                                                           weight=(1.0/(Grid_x+1.0))*(1.0/(Grid_y+1.0));
                                                           for(l=0;l<rown;l++){
                                                                               if(fabs(Zcoord-Zdata[l])<0.000001){
                                                                                                    rouz+=Wdata[l]*weight;
                                                                                                    }
                                                                               if(fabs(Zcoord-ZZdata[l])<0.000001){
                                                                                                     rouz_A+=WAdata[l]*weight;
                                                                                                     rouz_B+=WBdata[l]*weight;
                                                                                                     rouz_C+=WCdata[l]*weight;
                                                                                                     rouz_D+=WDdata[l]*weight;
                                                                                                     rouz_E+=WEdata[l]*weight;
                                                                                                     rouz_F+=WFdata[l]*weight;
                                                                                                     }
                                                                               }
                                                           fprintf(output3,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",Zcoord,rouz,rouz_A,rouz_B,rouz_C,rouz_D,rouz_E,rouz_F);
                                                           }

      fclose(output7);
      fclose(output3);
      fclose(output2);
      fclose(output4);
      fclose(output8); fclose(output9); fclose(output10);
      fclose(output11); fclose(output12);
      free(Xdata);
      free(Ydata);
      free(Zdata);
      free(Wdata);
      free(XXdata);
      free(YYdata);
      free(ZZdata);
      free(WAdata); free(WDdata);
      free(WBdata); free(WEdata);
      free(WCdata); free(WFdata);
      //Free the memory allocated to each array
      }
      free(super_x); free(super_y); free(super_z);
// Delete the transition file
      int status, Q;
      int status1, status2;
      char files2[128]; 
      char files3[128];
      char files4[128];
      printf("THE INTERMEDIATE FILES WILL BE REMOVED NOW\n");
      for(Q=1;Q<=frame;Q++){
                            sprintf(files2,"../out/MDS/Mass-%d.dat",Q);
                            sprintf(files3,"../out/MDS/Mass-%d-tmp.dat",Q);
                            sprintf(files4,"../out/MDS/Mass_species-%d-tmp.dat",Q);
                            status=remove(files2);
                            status1=remove(files3);
                            status2=remove(files4);
                            if(status==0&&status1==0&&status2==0){
                                          printf("%s file is deleted successfully.\n",files2);
                                          printf("%s file is deleted successfully.\n",files3);
                                          printf("%s file is deleted successfully.\n",files4);
                                          }
                            else {
                                  printf("Unable to delete the file\n");
                                  perror("ERROR");
                                  }
                            }
      
      printf("ALL CALCULATIONS ARE FINISHED\n");
      printf("RESULTS ARE SAVED IN THE OUTPUT FILES\n");                                                                                                                                                              
      return 0;
}

void error_g ()
{
	float g(double, double); // Mass smearing function
	int i;
	double r, g1, g2;
	printf ("r g(Fe) g(Cr)\n");
	for (i = 0; i < 10000; i ++) {
		r = i * 0.01;
		g1 = 55.84 * g(delta, r);
		g2 = 52. * g(delta, r);
		printf ("%f %e %e\n", r, g1, g2);
	}
	exit (1);
}

float f(double x, double y, double z)
{
      return sqrt(x*x+y*y+z*z);
      }           
float g(double x, double y)
{
      return pow(2.0*PI*x*x,-3.0/2)*exp(-y*y/(2.0*x*x));
      }             
