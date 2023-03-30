#ifndef ADJUSTMENT_ASSEMBLE_LINSYS_H_
#define ADJUSTMENT_ASSEMBLE_LINSYS_H_

#include "adjustment_impl.h"

/*
 * if vertex is known, return 1
 * */
int IsVertexKnown( int, int *, int );

void assemble_weight_matrix( double * *, double * *, int, int );
int lagrange_index_diag( int );
int lagrange_index_upper( int );
int lagrange_index_lower( int );

#endif
