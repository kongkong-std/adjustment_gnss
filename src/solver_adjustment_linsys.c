#include "../include/solver_adjustment_linsys.h"

void solver_adjust_linsys( double * * mat, double * rhs, int row, int column,
	double * sol )
{
    // solving least square equation
    /*
     * x = argmin || Ax - b ||
     * QR decomposition for A with Givens rotation
     * ( Ax - b ) ^ T ( Ax - b ) = () ^ T Q Q ^ T ()
     * Q ^ T A = R, Q = G1 G2 ... Gk
     * */
    givens_decomposition_impl( mat, rhs, row, column );

    // solution for least square equation
    /*
     * after QR decomposition
     * A is an upper triangular matrix
     * */
    solution_adjust_linsys( mat, rhs, column, sol );
}

void solution_adjust_linsys( double * * mat, double * rhs, int column, double * sol )
{
    // solving upper triangular linear system with direct method
    /*
     * x_k = ( b_k - \sum{ j = k + 1, n } a_kj x_j ) / a_kk
     * */
    for( int index = column - 1; index >= 0; index-- )
    {
	double sum = 0;
	for( int index_j = index + 1; index_j < column; index_j++ )
	{
	    sum += mat[ index ][ index_j ] * sol[ index_j ];
	}
	sol[ index ] = ( rhs[ index ] - sum ) / mat[ index ][ index ];
    }
}

void givens_decomposition_impl( double * * mat, double * rhs, int row, int column )
{
    for( int index = 0; index < column; index++ )
    {
	for( int index_j = index + 1; index_j < row; index_j++ )
	{
	    givens_rotation_impl( mat[ index ][ index ], mat[ index_j ][ index ], index, index_j,
		    row, column, mat, rhs );
	}
    }
}

void givens_rotation_impl( double val_a, double val_b, int index_i, int index_j,
	int row, int column, double * * mat, double * rhs )
{
    double givens_tau = 0, givens_c = 0, givens_s = 0;
    if( val_a != 0 && val_b != 0 )
    {
	if( fabs( val_b ) > fabs( val_a ) )
	{
	    givens_tau = - val_a / val_b;
	    givens_s = 1 / sqrt( 1 + givens_tau * givens_tau );
	    givens_c = givens_s * givens_tau;
	}
	else
	{
	    givens_tau = -val_b / val_a;
	    givens_c = 1 / sqrt( 1 + givens_tau * givens_tau );
	    givens_s = givens_c * givens_tau;
	}
    }
    else if( val_a != 0 && val_b == 0 )
    {
	givens_c = 1;
	givens_s = 0;
    }
    else if( val_a == 0 && val_b != 0 )
    {
	givens_c = 0;
	givens_s = 1;
    }
    else if( val_a == 0 && val_b == 0 )
    {
	givens_c = 1;
	givens_s = 0;
    }

    // computing mat
    /*
     * mat( index_i, index ) = c * mat( index_i, index ) - s * mat( index_j, index )
     * mat( index_j, index ) = s * mat( index_i, index ) + c * mat( index_j, index )
     * */
    double temp_1 = 0, temp_2 = 0;
    for( int index = 0; index < column; index++ )
    {
	temp_1 = mat[ index_i ][ index ];
	temp_2 = mat[ index_j ][ index ];
	mat[ index_i ][ index ] = givens_c * temp_1 - givens_s * temp_2;
	mat[ index_j ][ index ] = givens_s * temp_1 + givens_c * temp_2;
    }

    // computing rhs
    /*
     * rhs( index_i, 1 ) = c * rhs( index_i, 1 ) - s * rhs( index_j, 1 )
     * rhs( index_j, 1 ) = s * rhs( index_i, 1 ) + c * rhs( index_j, 1 )
     * */
    temp_1 = rhs[ index_i ];
    temp_2 = rhs[ index_j ];
    rhs[ index_i ] = givens_c * temp_1 - givens_s * temp_2;
    rhs[ index_j ] = givens_s * temp_1 + givens_c * temp_2;
}

void computing_residual( double * res, double * rhs, double * * mat, double * sol,
	int mat_row, int mat_column )
{
    // r = b - Ax
    /*
     * r_i = b_i - ( Ax )_i
     * ( Ax )_i = \sum _ j ^ n A_ij x_j
     * */
    for( int index = 0; index < mat_row; index++ )
    {
	double sum_temp = 0;
	for( int index_j = 0; index_j < mat_column; index_j++ )
	{
	    sum_temp += mat[ index ][ index_j ] * sol[ index_j ];
	}

	res[ index ] = rhs[ index ] - sum_temp;
    }
}
