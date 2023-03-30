#ifndef ADJUSTMENT_IMPL_H_
#define ADJUSTMENT_IMPL_H_

// header
#include "main.h"

/*
 * if initiali vertex and teminal vertex are known vertices,
 * return 1
 * */
int IsKnownVertex( int *, int *, int, int );

/*
 * data_mat, data_rhs, data_weight
 * data_number: baseline vertex number
 * data_src: baseline observation data
 * number_known: number of known vertex
 * data_known: coordinate of known vertex
 * data_row, count_known: known vertex
 * count: count of vertex
 * row_count_unknown
 * */
void assemble_adjust_linsys( double * *, double *, double * *,
	int * *, double * *,
	int *, double * *,
	int, int, int, int );

void update_linsys_weight_mat( double * *, double * *, double *, int, int );

void direct_solver_linsys( double * *, double *, double *, int );
void thomas_algorithm_linsys( double *, double *, double *, double *, double *, int );

void solver_adjust_linsys( double * *, double *, int, int, double * );

#endif
