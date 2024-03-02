## Overview

ScaleLat (Scale Lattice) is a computer program written in C ++ for performing the atomic structure analysis of multi-phase system or high entropy alloys (HEAs). The program implements an atomic cluster extraction algorithm to obtain all independent and symmetry-reduced characteristic chemical structures for the complex atomic configurations which are usually obtained from molecular dynamics or kinetic Monte-Carlo simulations for supercell containing more than 104 atoms. ScaleLat employes an efficient and unique chemical structure matching algorithm to map all extracted atomic clusters from a large supercell (>104 atoms) to a representative small one (~ 103 or less), providing the possibility to directly use the highly accurate quantum mechanic methods to study the electronic, magnetic, and mechanical properties of multi-component alloys with complex microstructures.

Nature of Problem: Very large supercells containing more than 104 atoms must be used to describe the atomic structure of complex multi-phase or high entropy alloys, and their properties are almost not possible to investigate using quantum mechanical calculations within the conventional computational facilities right now. Decreasing the total number of atoms in a smaller representing supercell which preserves all essential chemical structures of the original large supercell is the key to the problem. 

Solution to the Problem: An atomic cluster extraction algorithm is realized to thoroughly analysis all essential chemical structures of multi-phase system or high entropy alloys; Utilizing the direct atom swapping method, the chemical structure matching procedure is developed to map all extracted characteristic atomic clusters of benchmark structure to the small supercell by minimizing their differences in both the types and relative proportions of atomic cluster sets. 


## software environment
Only supported by Ubuntu 16. The g++, gfortran and cmake compilers need to be installed in advance.

    $ sudo apt-get update
    $ sudo apt-get install g++
    $ sudo apt-get install gfortran
    $ sudo apt-get install cmake
		
The software OVITO can be used to visualize the results files "*.movie" or "*.xyz". 

    $ sudo apt-get install ovito

## installation
Download the file ScaleLat.zip.

## Run the following commands to extract and access to all of the original files:
    $ unzip ScaleLat.zip 
    $ sudo chmod -R 775 ScaleLat
    $ cd ScaleLat

## Configure, compile, and install ScaleLat with:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

## How to run examples:

    $ cd example/the/folder/contains/*.sh
    $ bash *.sh

## Run FEcMD

There are various commands to perform the corresponding function.
Perform kMC simulation:

    $ ./scalelat kMC

Perform cluster extraction:

    $ ./scalelat CE

Perform chemical structure matching:

    $ ./scalelat MAP

Perform cluster analysis for density, size and volume of precipitates:

    $ ./scalelat analysis

Perform coarse-grained analysis to get mass density distribution:

    $ ./scalelat MDS

Calculate SRO

    $ ./scalelat SRO

Calculate ACF

    $ ./scalelat ACF

## input file:
-------------parameters for run.in (## ... should be deleted in run.in):
eam_file	../in/potential/Fe_Cr_Eich_2015_TBM_lammps.eam.alloy  ## the atomic potential funtion
MC_method	vacancy			## vacancy or atom exchange
cellSize	10 10 10		## the size of cell expansion
cell_order	order  			## order or disorder after cell expansion
threshold	0.001			## threshold, default: 0.001
NeighbourMax	20			## the max neighbour number for the xth shell
T_max		673			## temperature [K]
Nelems		2			## the number of elements
elemName	Fe Cr			## elements names
NLOOP		2.e10			## loop number
stepAvg		1.e5			## average steps for calculating energy 
stepMovie	1.e7			## average steps for output *.movie
stepVacancy	100			## the steps of the vacancy exchange randomly, default: 100
stepAtom	1			## the steps of the atom exchange randomly, default: 1
nthreads	16			## the number of cores that run in parallel
seed		7			## random seed
nebrTabFac	200			## the max neighbour number within cutoff radius
--------------------------------------------------------------
----------------------------RDF-------------------------------
limitRdf	100			# The limit steps of RDF
stepRdf		100000			# The average steps of RDF
rangeRdf	5.			# The range of RDF [Angstrom]
sizeHistRdf	200			# The size of histogram
--------------------------------------------------------------
-----------------------------MAP------------------------------
clusterSize	6 6 6			# The size of cluster extraction [unit: lattice parameter]
error_map	1			# the standard deviation (1) or the mean absolute error (0)
pointGroup	m3m			# the symmetry points group
cell_order	disorder		# order or disorder for the initial trial supercell
--------------------------------------------------------------
--------------------------Analysis----------------------------
ith_nebrR	2			# The number of shells of cluster analysis
C_threshold	0.95			# The concentration threshold of cluster analysis
elem_preci	Cr			# The element type of cluster analysis

Alloy atomic trajectory input file ‘benchmark.xyz’ format (example for HEA):
108000  ## atomic number
Lattice="105.600002 0.0 0.0 0.0 105.599999 0.0 0.0 0.0 105.599999" Properties=species:S:1:pos:R:3  ## box region and the title for ovito
Cr 44.0 10.56 12.320001
Cr 82.720001 95.040001 79.199997
Fe 77.440002 33.439999 15.84
Co 22.880001 88.0 29.92
Co 40.48 26.4 59.84
Ni 44.0 0.0 82.720001
… ## The following lines are 4 arrays of tabulated values: element_type position_x position_y position_z

The unit cell ‘CONTCAR’ format (example for HEA):
POSCAR file written by OVITO  ## title
1.0  ## scaling factor
       26.3999996185         0.0000000000         0.0000000000
        0.0000000000        28.1599998474         0.0000000000
        0.0000000000         0.0000000000        28.1599998474  ## box region
   Co   Cr   Cu   Fe   Ni  ## elements types
  364  402  435  368  351  ## atomic number
Direct  ## must be fractional coordinates
     0.200000003         0.000000000         0.812500000
     0.333333343         0.437500030         0.000000000
     0.866666734         0.125000000         0.562500000
     0.733333349         0.500000000         0.937500000
     0.066666663         0.250000000         0.562500000
     0.733333349         0.750000060         0.187500015
…  ## The following lines are 3 arrays of tabulated values: position_x position_y position_z

## output file:

*.xyz and *.movie files can be visualized by ‘ovito’.
md: this folder includes atomic trajectory files.
supercell: this folder includes cell expansion files.
sro.dat: SRO data.
rdf.dat: RDF data.
precipitates_*.dat: the average size, density and volume of precipitates.
precipitates_*.movie: the atom structure of precipitates.
CE_benchmark: this folder includes extracted clusters of the benchmark supercell.
CE_trial: this folder includes extracted clusters of the small trial structure.
CE/CE.dat: the number of characteristic chemical clusters.
CE/CE*.xyz: characteristic chemical clusters structure.
CE/concentration.dat: the atomic concentration of characteristic chemical clusters.
diff.dat: STD (or MAE) of atomic concentration and clusters proportion pi.
mapping.dat: characteristic chemical clusters proportion pi of the benchmark supercell and the small trial supercell.
mapping.xyz: atomic structure file of the matching structure.
