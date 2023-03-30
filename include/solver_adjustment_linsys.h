#ifndef SOLVER_ADJUSTMENT_LINSYS_H_
#define SOLVER_ADJUSTMENT_LINSYS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void givens_decomposition_impl( double * *, double *, int, int );
void givens_rotation_impl( double, double, int, int, int, int,
	double * *, double * );
void solution_adjust_linsys( double * *, double *, int, double * );

#endif
