#ifndef ADJUSTMENT_IMPL_H_
#define ADJUSTMENT_IMPL_H_

#include <stdio.h>
#include <stdlib.h>

void assemble_adjust_linsys( int * *, double * *, int, int,
	double *, double *, int, int,
	double * *, double *, double * *);

void update_linsys_weight_mat( double * *, double * *, double *, int, int );

void direct_solver_linsys( double * *, double *, double *, int );
void thomas_algorithm_linsys( double *, double *, double *, double *, double *, int );

void solver_adjust_linsys( double * *, double *, int, int, double * );

#endif
