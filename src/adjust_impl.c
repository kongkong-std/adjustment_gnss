#include "adjust_impl.h"

void adjustment_implementation( double * * src_data, int size_row, int size_column )
{
    // solution
    /*
     * sol_data coordinate value
     * id 1, id 2 known value
     * id 3 coordinate x y z
     * id 4 coordinate x y z
     * id 5 coordinate x y z
     * id 6 coordinate x y z
     * */
    double * * sol_data = NULL;
    if( ( sol_data = ( double * * ) malloc( 4 * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < 4; index++ )
    {
	if( ( *( sol_data + index ) = ( double * ) malloc( 3 * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // free memory
    for( int index = 0; index < 4; index++ )
    {
	free( *( sol_data + index ) );
    }
    free( sol_data );
}
