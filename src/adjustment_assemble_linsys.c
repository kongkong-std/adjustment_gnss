#include "../include/adjustment_assemble_linsys.h"

void assemble_adjust_linsys( double * * mat, double * rhs, double * * weight_mat,
	int * * data_number, double * * data_src,
	int * number_known, double * * data_known,
	int data_row, int count_known,
	int count, int row_count_unknown,
	double * * original_mat, double * original_rhs )
{
    puts( "============adjustment linear system assembling============" );

    // cases
    /*
     * baseline
     * unknown node -- unknown node
     * known node -- known node
     * unknown node -- known node
     * */

    // assemble mat, weight_mat, rhs
    for( int index = 0; index < data_row; index++ )
    {
	if( IsVertexKnown( data_number[ index ][ 0 ], number_known, count_known ) != 1 
		&& IsVertexKnown( data_number[ index ][ 1 ], number_known, count_known ) != 1 )
	{
	    // unknown node -- unknown node
	    for( int index_i = 0; index_i < 3; index_i++ )
	    {
		weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index_diag( index_i ) ];
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] = 1;
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] = -1;
		rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ];

		// original linear system
		original_mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] = 1;
		original_mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] = -1;
		original_rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ];

		// mat = weight_mat * mat, rhs = weight_mat * rhs
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] *= 
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] *= 
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
	    }
	}
	else if( IsVertexKnown( data_number[ index ][ 0 ], number_known, count_known ) == 1 
		&& IsVertexKnown( data_number[ index ][ 1 ], number_known, count_known ) == 1 )
	{
	    // known node -- known node
	    // keep mat, rhs
	}
	else
	{
	    // known node -- unknown node
	    if( IsVertexKnown( data_number[ index ][ 0 ], number_known, count_known ) != 1 )
	    {
		// data_number[ index ][ 0 ] is unknown node, data_number[ index ][ 1 ] is known node
		for( int index_i = 0; index_i < 3; index_i++ )
		{
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index_diag( index_i ) ];
		    mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] = -1;
		    rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] 
			- data_known[ data_number[ index ][ 1 ] ][ index_i ];

		    // original linear system
		    original_mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] = -1;
		    original_rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] 
			- data_known[ data_number[ index ][ 1 ] ][ index_i ];

		    // mat = weight_mat * mat, rhs = weight_mat * rhs
		    mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - count_known ) + index_i ] *=
		       	weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		    rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		}
	    }
	    else
	    {
		// data_number[ index ][ 1 ] is unknown node, data_number[ index ][ 0 ] is known node
		for( int index_i = 0; index_i < 3; index_i++ )
		{
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index_diag( index_i ) ];
		    mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] = 1;
		    rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] 
			+ data_known[ data_number[ index ][ 0 ] ][ index_i ];

		    // original linear system
		    original_mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] = 1;
		    original_rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] 
			+ data_known[ data_number[ index ][ 0 ] ][ index_i ];

		    // mat = weight_mat * mat, rhs = weight_mat * rhs
		    mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - count_known ) + index_i ] *=
		       	weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		    rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		}
	    }
	}
    }
}

int IsVertexKnown( int value, int * number_known, int size_number_known )
{
    int status = 0;

    for( int index = 0; index < size_number_known; index++ )
    {
	if( *( number_known + index ) == value )
	{
	    status = 1;
	    break;
	}
    }

    return status;
}

void assemble_weight_matrix( double * * weight_mat, double * * data_src, int index, int index_i )
{
		// block weight_mat assembling
		/*
		 * upper triangular matrix, index_i = 0, 1, 2, index_j >= index_i
		 * weight_mat[ 3 * index + index_i ][ 3 * index + index_j ] = data_src[ index ][ 3 * index_i + index_j  ]
		 *
		 * lower triangular matrix, index_i = 1, 2, index_j < index_i
		 * symmetric matrix
		 * */
		// diagonal element
		weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = data_src[ index ][ lagrange_index_diag( index_i ) ];

		// upper triangular matrix
		for( int index_j = index_i + 1; index_j < 3; index_j++ )
		{
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_j ] = data_src[ index ][ lagrange_index_upper( 3 * index_i + index_j ) ];
		}

		// lower triangular matrix
		for( int index_j = 0; index_j < index_i; index_j++ )
		{
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_j ] = data_src[ index ][ lagrange_index_lower( 3 * index_i + index_j ) ];
		}

}

int lagrange_index_lower( int x )
{
    int value = 0;

    // f(3) = 4, f(6) = 5, f(7) = 7
    int f_3 = 4, f_6 = 5, f_7 = 7;

    value = f_3 * ( x - 6 ) * ( x - 7 ) / ( 3 - 6 ) / ( 3 - 7 )
	+ f_6 * ( x - 3 ) * ( x - 7 ) / ( 6 - 3 ) / ( 6 - 7 )
	+ f_7 * ( x - 3 ) * ( x - 6 ) / ( 7 - 3 ) / ( 7 - 6 );

    return value;
}

int lagrange_index_upper( int x )
{
    int value = 0;

    // f(1) = 4, f(2) = 5, f(5) = 7
    int f_1 = 4, f_2 = 5, f_5 = 7;

    value = f_1 * ( x - 2 ) * ( x - 5 ) / ( 1 - 2 ) / ( 1 - 5 )
	+ f_2 * ( x - 1 ) * ( x - 5 ) / ( 2 - 1 ) / ( 2 - 5 )
	+ f_5 * ( x - 1 ) * ( x - 2 ) / ( 5 - 1 ) / ( 5 - 2 );

    return value;
}

int lagrange_index_diag( int x )
{
    int value = 0;
    int f_0 = 3, f_1 = 6, f_2 = 8;

    value = f_0 * ( x - 1 ) * ( x - 2 ) / 2
	+ f_1 * x * ( x - 2 ) / ( -1 )
	+ f_2 * x * ( x - 1 ) / 2;

    return value;
}

void update_linsys_weight_mat( double * * weight_mat, double * * mat, double * rhs, int row, int column )
{
    // updating linsys
    /*
     * mat = inverse( weight_mat ) * mat
     * rhs = inverse( weight_mat ) * rhs
     *
     * size weight_mat = row x row
     * size mat = row x column
     * size rhs = row x 1
     * */

    // size pointer_temp = row x 1
    double * pointer_temp = NULL;
    if( ( pointer_temp = ( double * ) malloc( row * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize pointer_temp = 0
    for( int index = 0; index < row; index++ )
    {
	*( pointer_temp + index ) = 0;
    }

    // updating rhs
    /*
     * solving linear system
     * weight_mat * pointer_temp = rhs
     * */
    direct_solver_linsys( weight_mat, rhs, pointer_temp, row );
    for( int index = 0; index < row; index++ )
    {
	*( rhs + index ) = *( pointer_temp + index );
    }

    // updating mat
    /*
     * solving linear system
     * weight_mat * pointer_temp = mat( :, j ), j = 1, 2, ..., column
     * */
    double * mat_column = NULL;    // size mat_column = row x 1
    if( ( mat_column = ( double * ) malloc( row * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    for( int index = 0; index < column; index++ )
    {
	for( int index_i = 0; index_i < row; index_i++ )
	{
	    *( mat_column + index_i ) = mat[ index_i ][ index ];
	}

	direct_solver_linsys( weight_mat, mat_column, pointer_temp, row );

	for( int index_i = 0; index_i < row; index_i++ )
	{
	    mat[ index_i ][ index ] = *( pointer_temp + index_i );
	}
    }

    // free
    free( mat_column );
    free( pointer_temp );
}

void direct_solver_linsys( double * * A, double * b, double * x, int size )
{
    // direct method solving linear system
    /*
     * A x = b, size x size
     *
     * A is tri-diagonal matrix, use Thomas algorithm
     * */
    // diagonal, up sub-diagonal, low sub-diagonal
    /*
     * size lower = size - 1 x 1
     * size diag = size x 1
     * size upper = size -1 x 1
     * */
    double * lower = NULL, * diag = NULL, * upper = NULL;
    if( ( lower = ( double * ) malloc( ( size - 1 ) * sizeof( double ) ) ) == NULL 
	    || ( diag = ( double * ) malloc( size * sizeof( double ) ) ) == NULL 
	    || ( upper = ( double * ) malloc( ( size - 1 ) * sizeof( double ) ) ) == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize lower = 0, diag = 0, upper = 0
    for( int index = 0; index < size - 1; index++ )
    {
	*( lower + index ) = 0;
	*( diag + index ) = 0;
	*( upper + index ) = 0;
    }
    *( diag + size - 1 ) = 0;

    // set value
    /*
     * lower = A( i + 1, i ), i = 1, ..., size - 1
     * diag = A( i, i ), i = 1, 2, ..., size
     * upper = A( i, i + 1 ), i = 1, 2, ..., size -1
     * */
    for( int index = 0; index < size - 1; index++ )
    {
	*( lower + index ) = A[ index + 1 ][ index ];
	*( diag + index ) = A[ index ][ index ];
	*( upper + index ) = A[ index ][ index + 1 ];
    }
    *( diag + size - 1 ) = A[ size - 1 ][ size - 1 ];

    thomas_algorithm_linsys( lower, diag, upper, b, x, size );

    // free
    free( lower );
    free( diag );
    free( upper );
}

void thomas_algorithm_linsys( double * a, double * d, double * c, double * b, double * x, int n )
{
    // parameter
    /*
     * lower = [ a_1, a_2, ..., a_{n-1} ]
     * diag = [ d_1, d_2, ..., d_n ]
     * upper = [ c_1, c_2, ..., c_{n-1} ]
     * rhs = [ b_1, b_2, ..., b_{n-1} ]
     * solution = [ x_1, x_2, ..., x_{n-1} ]
     * */

    // temporary variable
    /*
     * size temp_p = n x 1
     * size temp_q = n - 1 x 1
     * size temp_y = n x 1
     * */
    double * temp_p = NULL, * temp_q = NULL, * temp_y = NULL;
    if( ( temp_p = ( double * ) malloc( n * sizeof( double ) ) ) == NULL ||
	    ( temp_q = ( double * ) malloc( ( n - 1 ) * sizeof( double ) ) ) == NULL ||
	    ( temp_y = ( double * ) malloc( n * sizeof( double ) ) ) == NULL)
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize temp_p = 0, temp_q = 0, temp_y = 0
    for( int index = 0; index < n - 1; index++ )
    {
	*( temp_p + index ) = 0;
	*( temp_q + index ) = 0;
	*( temp_y + index ) = 0;
    }
    *( temp_p + n - 1 ) = 0;
    *( temp_y + n - 1 ) = 0;

    // computing temp_p, temp_q
    /*
     * p_1 = d_1, q_1 = c_1 / d_1
     * p_k = d_k - a_k q_{k-1}, q_k = c_k / q_k, k = 2, ..., n - 1
     * p_n = d_n - a_n q_{n-1}
     * */
    *temp_p = *d;
    *temp_q = *c / *d;
    for( int index = 1; index < n - 1; index++ )
    {
	*( temp_p + index ) = *( d + index ) - *( a + index - 1 ) 
	    * ( *( temp_q + index - 1 ) );
	* ( temp_q + index ) = *( c + index ) / ( *( temp_p + index ) );
    }
    *( temp_p + n - 1 ) = *( d + n - 1 ) - *( a + n - 1 ) * ( *( temp_q + n - 2 ) );

    // computing temporary solution
    /*
     * y_1 = b_1 / d_1
     * y_k = ( b_k - a_k y_{k-1} ) / p_k, k = 2, 3, ..., n
     * */
    *temp_y = *b / *d;
    for( int index = 1; index < n; index++ )
    {
	*( temp_y + index ) = ( *( b + index ) - *( a + index ) * ( *( temp_y + index - 1 ) ) ) / ( *( temp_p + index ) );
    }

    // computing solution
    /*
     * x_n = y_n
     * x_k = y_k - q_k x_{k+1}, k = n - 1, n - 2, ..., 1
     * */
    *( x + n - 1 ) = *( temp_y + n - 1 );
    for( int index = n - 2; index >= 0; index-- )
    {
	*( x + index ) = *( temp_y + index ) - *( temp_q + index ) * ( *( x + index + 1 ) );
    }

    // free
    free( temp_p );
    free( temp_q );
    free( temp_y );
}
