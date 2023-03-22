#ifndef ADJUSTMENT_IMPL_H_
#define ADJUSTMENT_IMPL_H_

#include <stdio.h>
#include <stdlib.h>

void assemble_adjust_linsys( int * *, double * *, int, int,
	double *, double *, int, int,
	double * *, double *, double * *);
void solver_adjust_linsys( double * *, double *, int, int, double * );

#endif
