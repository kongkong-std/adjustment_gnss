#include "adjustment_impl.h"

void adjustment_impl( int * * data_number, double * * data_src,
	int data_row, int data_column, int data_known )
{
    puts( "============adjustment of GNSS network implementation============" );

    /*
     * number of known data
     *     node number 1
     *     node number 2
     * */
    int known_number_1 = 1, known_number_2 = 2;
    
    //double known_solution_1[ 3 ] = { 402.35087, -4652995.30109, 4349760.77753 };
    //double known_solution_2[ 3 ] = { 8086.03178, -4642712.84739, 4360439.08326 };
    
    double known_solution_1[ 3 ] = { 3213503.9707, 852378.2675, 402.9194 }; 
    double known_solution_2[ 3 ] = { 3213519.0995, 852411.9884, 387.2503 }; 

    /*
     * count of observation data contains unknown nodes
     * */
    int count_unknown = 0;
    for( int index_i = 0; index_i < data_row; index_i++ )
    {
#if 1
	if( ! ( ( ( data_number[ index_i ][ 0 ] == known_number_1 ) && ( data_number[ index_i ][ 1 ] == known_number_2 ) ) 
		|| ( ( data_number[ index_i ][ 0 ] == known_number_2 ) && ( data_number[ index_i ][ 1 ] == known_number_1 ) ) ) )
	{
	    count_unknown++;
	}
#endif

#if 0
	if( !( data_number[ index_i ][ 0 ] == known_number_1 && data_number[ index_i ][ 1 ] == known_number_1 ) )
	{
	    count_unknown++;
	}
#endif
    }
#if 1
    printf( ">>>>count of unknown is %d\n", count_unknown );
#endif

    /*
     * array of unknown node
     * */
    int * node_unknown = NULL;
    int count_node_unknown = data_row * 2;
    if( ( node_unknown = ( int * ) malloc( count_node_unknown * sizeof( int ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < 2; index++ )
    {
	for( int index_i = 0; index_i < data_row; index_i++ )
	{
	    node_unknown[ index * data_row + index_i ] = data_number[ index_i ][ index ];
	}
    }
    for( int index_i = 0; index_i < count_node_unknown; index_i++ )
    {
	for( int index_j = index_i + 1; index_j < count_node_unknown; index_j++ )
	{
	    if( node_unknown[ index_j ] == node_unknown[ index_i ] )
	    {
		// deleting node_unknown[ index_j ]
		for( int index_k = index_j; index_k < count_node_unknown; index_k++ )
		{
		    node_unknown[ index_k ] = node_unknown[ index_k + 1 ];
		}
		count_node_unknown--;
		index_j--;
	    }
	}
    }
    for( int index = 0; index < count_node_unknown; index++ )
    {
	if( node_unknown[ index ] == known_number_1 || node_unknown[ index ] == known_number_2 )
	//if( node_unknown[ index ] == known_number_1 )
	{
	    // delete known node
	    for( int index_i = index; index_i < count_node_unknown; index_i++ )
	    {
		node_unknown[ index_i ] = node_unknown[ index_i + 1 ];
	    }
	    count_node_unknown--;
	    index--;
	}
    }
#if 1    // print unknown node
    printf( ">>>>number of unknown node is %d\n", count_node_unknown );
    for( int index = 0; index < count_node_unknown; index++ )
    {
	printf( "%d\n", node_unknown[ index ] );
    }
#endif

    // sorting node in number order
    for( int index = 0; index < count_node_unknown; index++ )
    {
	for( int index_i = index + 1; index_i < count_node_unknown; index_i++ )
	{
	    if( *( node_unknown + index_i ) < *( node_unknown + index ) )
	    {
		int temp = 0;
		temp = *( node_unknown + index_i );
		*( node_unknown + index_i ) = *( node_unknown + index );
		*( node_unknown + index ) = temp;
	    }
	}
    }
#if 1    // print sorted node
    printf( ">>>>after sorting...\n" );
    for( int index = 0; index < count_node_unknown; index++ )
    {
	printf( "%d\n", *( node_unknown + index ) );
    }
#endif

    // linear system
    /*
     * P_j = P_i + P_ij, i \ne j
     * P_i = [x_i, y_i, z_i], i = 1, 2, ..., node in GNSS network
     * P_ij observation data
     *
     * A x = b
     *     size A = ( count_unknown * 3 ) x ( count_node_unknown * 3 )
     *
     * right-hand side b
     *     size b = ( count_unknown * 3 ) x 1
     *
     * solution x
     *     size x = ( count_node_unknown * 3 ) x 1
     *
     * weighted matrix W_i
     *     W_i is diagonal matrix, size W_i = 3 x 3
     *     W_i = diag( sigma_xx, sigma_yy, sigma_zz );
     *     W is block diagonal matrix, W_1, W_2, ..., size W = ( count_unknown * 3 ) x ( count_unknown * 3 )
     * */
    int data_mat_row = count_unknown * 3, data_mat_column = count_node_unknown * 3;
    double * * data_mat = NULL;
    double * * data_weight = NULL;
    double * data_rhs = NULL;
    double * data_sol = NULL;

    // data_mat
    /*
     * size data_mat = data_mat_row x data_mat_column
     * */
    if( ( data_mat = ( double * * ) malloc( data_mat_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < data_mat_row; index++ )
    {
	if( ( *( data_mat + index ) = ( double * ) malloc( data_mat_column * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // data_weight
    /*
     * size data_weight = data_mat_row x data_mat_row
     * */
    if( ( data_weight = ( double * * ) malloc( data_mat_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < data_mat_row; index++ )
    {
	if( ( *( data_weight + index ) = ( double * ) malloc( data_mat_row * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // data_rhs
    /*
     * size data_rhs = data_mat_row x 1
     * */
    if( ( data_rhs = ( double * ) malloc( data_mat_row * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // data_sol
    /*
     * size data_sol = data_mat_column x 1
     * */
    if( ( data_sol = ( double * ) malloc( data_mat_column * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize data_mat = 0
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_column; index_j++ )
	{
	    data_mat[ index_i ][ index_j ] = 0;
	}
    }

#if 1    // print data_mat
    printf( ">>>>mat_row = %d\n>>>>mat_column = %d\n", data_mat_row, data_mat_column );
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_column; index_j++ )
	{
	    printf( "%lf ", data_mat[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif

    // initialize data_weight = 0
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_row; index_j++ )
	{
	    data_weight[ index_i ][ index_j ] = 0;
	}
    }
#if 1    // print data_weight
    printf( ">>>>weight_row = %d\n>>>>weight_column = %d\n", data_mat_row, data_mat_row );
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_row; index_j++ )
	{
	    printf( "%lf ", data_weight[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif

    // initialize data_rhs = 0
    for( int index = 0; index < data_mat_row; index++ )
    {
	data_rhs[ index ] = 0;
    }

#if 1    // print data_rhs
    for( int index = 0; index < data_mat_row; index++ )
    {
	printf( "%lf\n", data_rhs[ index ] );
    }
#endif

    // initialize data_sol = 0
    for( int index = 0; index < data_mat_column; index++ )
    {
	data_sol[ index ] = 0;
    }

#if 0    // print data_sol
    for( int index = 0; index < data_mat_column; index++ )
    {
	printf( "%lf\n", data_sol[ index ] );
    }
#endif

    // assemble matrix and rhs
    /*
     * known node -- unknown node
     * unknown node -- known node
     * */
    assemble_adjust_linsys( data_number, data_src, data_row, data_column,
	    known_solution_1, known_solution_2,
	    known_number_1, known_number_2,
	    data_mat, data_rhs, data_weight );
#if 1    // print data_mat
    printf( ">>>>mat_row = %d\n>>>>mat_column = %d\n", data_mat_row, data_mat_column );
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_column; index_j++ )
	{
	    printf( "%lf ", data_mat[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif
#if 1    // print data_rhs
    printf( ">>>>size data_rhs = %d\n", data_mat_row );
    for( int index = 0; index < data_mat_row; index++ )
    {
	printf( "%lf\n", data_rhs[ index ] );
    }
#endif
#if 1    // print data_weight
    printf( ">>>>weight_mat row = %d\n>>>>weight_mat column = %d\n", data_mat_row, data_mat_row );
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_row; index_j++ )
	{
	    printf( "%.2le ", data_weight[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif

    // updating linear system with weighted matrix
    /*
     * x = argmin || W A x - W b ||
     * W = inverse( data_weight )
     * */
   // update_linsys_weight_mat( data_weight, data_mat, data_rhs, data_mat_row, data_mat_column );

    // solver adjustment linear system
    /**/
    solver_adjust_linsys( data_mat, data_rhs, data_mat_row, data_mat_column,
	    data_sol );
#if 1    // print data_mat
    printf( ">>>>mat_row = %d\n>>>>mat_column = %d\n", data_mat_row, data_mat_column );
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_column; index_j++ )
	{
	    printf( "%lf ", data_mat[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif
#if 1    // print data_rhs
    printf( ">>>>size data_rhs = %d\n", data_mat_row );
    for( int index = 0; index < data_mat_row; index++ )
    {
	printf( "%lf\n", data_rhs[ index ] );
    }
#endif
#if 1    // print data_sol
    printf( ">>>>size data_sol = %d\n", data_mat_column );
    for( int index = 0; index < data_mat_column; index++ )
    {
	printf( "%lf\n", data_sol[ index ] );
    }
#endif

    // free
    free( data_sol );
    free( data_rhs );
    for( int index = 0; index < data_mat_row; index++ )
    {
	free( *( data_mat + index ) );
    }
    free( data_mat );
    free( node_unknown );
}
