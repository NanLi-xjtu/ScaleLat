#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>

#include "define.h"
#include "print.h"
#include "superCell.h"
#include "mc.h"
#include "eamforce.h"
#include "eam_alloy.h"
#include "nebrList.h"
#include "measurement.h"
#include "Mass_density_distribution.h"
#include "mapping.h"
#include "symmetry.h"

extern char mode[128];
void Run (int argc, char **argv);

#endif
