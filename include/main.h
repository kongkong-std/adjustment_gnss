#ifndef MAIN_H_
#define MAIN_H_

// header
#include <stdio.h>
#include <stdlib.h>
#include "graph_impl.h"

// function prototype
/*
 * graph
 * count of known vertex
 * number of known vertex array
 * coordinate of known vertex array
 * */
void graph_adjustment( AdjGraph *, int, int *, double * * );

/*
 * adjustment implementation
 *     number of vertex
 *     observation data: baseline solution, co-variance
 *     number of known vertex
 *     position of known vertex
 *     row of observation data
 *     column of observation data
 *     count of known vertex
 *     count of vertex
 * */
void adjustment_impl( int * *, double * *, int *, double * *,
	int, int, int, int );

#endif
