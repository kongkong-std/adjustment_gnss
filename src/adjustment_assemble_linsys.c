#include "adjustment_assemble_linsys.h"

void assemble_adjust_linsys( int * * data_number, double * * data_src, int data_row, int data_column,
	double * known_solution_1, double * known_solution_2,
	int known_number_1, int known_number_2,
	double * * mat, double * rhs)
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
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
		mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
		rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ];
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
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] - known_solution_1[ index_i ];
		    }
		}
		else
		{
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 0 ] - 3 ) + index_i ] = -1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] - known_solution_2[ index_i ];
		    }
		}
	    }
	    if( data_number[ index ][ 1 ] != known_number_1 && data_number[ index ][ 1 ] != known_number_2 )
	    {
		// data_number[ index ][ 1 ] is unknown node, data_number[ index ][ 0 ] is known node
		if( data_number[ index ][ 0 ] == known_number_1 )
		{
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] + known_solution_1[ index_i ];
		    }
		}
		else
		{
		    for( int index_i = 0; index_i < 3; index_i++ )
		    {
			mat[ 3 * index + index_i ][ 3 * ( data_number[ index ][ 1 ] - 3 ) + index_i ] = 1;
			rhs[ 3 * index + index_i ] = data_src[ index ][ index_i ] + known_solution_2[ index_i ];
		    }
		}
	    }
	}
    }
}
