#include "adjustment_assemble_linsys.h"

void assemble_adjust_linsys( int * * data_number, double * * data_src, int data_row, int data_column,
	double * known_solution_1, double * known_solution_2,
	int known_number_1, int known_number_2,
	double * * mat, double * rhs, double * * weight_mat )
{
    puts( "============adjustment linear system assembling============" );

    // cases
    /*
     * unknown node -- unknown node
     * known node -- known node
     * known node -- unknown node
     * */
    for( int index = 0; index < data_row; index++ )
    {
	if( data_number[ index ][ 0 ] != known_number_1 && data_number[ index ][ 0 ] != known_number_2
	 && data_number[ index ][ 1 ] != known_number_1 && data_number[ index ][ 1 ] != known_number_2 )
	{
	    // unknown node -- unknown node
	    for( int index_i = 0; index_i < 3; index_i++ )
	    {
		weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index( index_i ) ];
		
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
		rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ];

#if 1
		// mat = weight_mat * mat, rhs = weight_mat * rhs
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] *= 
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] *= 
		    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
		rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
#endif
	    }
	}
	else if( ( data_number[ index ][ 0 ] == known_number_1 && data_number[ index ][ 1 ] == known_number_2 ) 
		|| ( data_number[ index ][ 0 ] == known_number_2 && data_number[ index ][ 1 ] == known_number_1 ) )
	{
	    // known node -- known node
	}
	else
	{
	    // known node -- unknown node
	    if( data_number[ index ][ 0 ] != known_number_1 && data_number[ index ][ 0 ] != known_number_2 )
	    {
		// data_number[ index ][ 0 ] is unknown node, data_number[ index ][ 1 ] is known node
		if( data_number[ index ][ 1 ] == known_number_1 )
		{
		    // known node is known_number_1
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
		        weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index( index_i ) ];
			
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] - known_solution_1[ index_i ];

#if 1
			// mat = weight_mat * mat, rhs = weight_mat * rhs
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] *=
			    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
			rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
#endif
		    }
		}
		else
		{
		    // known node is known_number_2
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
		        weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index( index_i ) ];
			
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] - known_solution_2[ index_i ];

#if 1
			// mat = weight_mat * mat, rhs = weight_mat * rhs
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] *=
			    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
			rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
#endif
		    }
		}
	    }
	    else if( data_number[ index ][ 1 ] != known_number_1 && data_number[ index ][ 1 ] != known_number_2 )
	    {
		// data_number[ index ][ 1 ] is unknown node, data_number[ index ][ 0 ] is known node
		if( data_number[ index ][ 0 ] == known_number_1 )
		{
		    // known node is known_number_1
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
		        weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index( index_i ) ];
			
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] + known_solution_1[ index_i ];

#if 1
			// mat = weight_mat * mat, rhs = weight_mat * rhs
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] *=
			    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
			rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
#endif
		    }
		}
		else
		{
		    // known node is known_number_2
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
		        weight_mat[ 3 * index + index_i ][ 3 * index + index_i ] = 1. / data_src[ index ][ lagrange_index( index_i ) ];
			
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] + known_solution_2[ index_i ];

#if 1
			// mat = weight_mat * mat, rhs = weight_mat * rhs
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] *=
			    weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
			rhs[ 3 * index + index_i ] *= weight_mat[ 3 * index + index_i ][ 3 * index + index_i ];
#endif
		    }
		}
	    }
	}
    }
}

int lagrange_index( int x )
{
    int value = 0;
    int f_0 = 3, f_1 = 6, f_2 = 8;

    value = f_0 * ( x - 1 ) * ( x - 2 ) / 2
	+ f_1 * x * ( x - 2 ) / ( -1 )
	+ f_2 * x * ( x - 1 ) / 2;

    return value;
}
