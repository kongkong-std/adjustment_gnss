#include "../include/adjustment_impl.h"

void graph_adjustment( AdjGraph * graph, int count_known, int * number_known, double * * data_known )
{
    int data_row = 0, data_column = 0;
    data_row = graph->count_edge;
    data_column = sizeof( graph->array->head->weight ) / sizeof( graph->array->head->weight[ 0 ] );

    int * * data_number = NULL;    // vertex number
    double * * data_src = NULL;    // observation data

    /*
     * size data_number = data_row x 2
     * size data_number = data_row x data_column
     * */
    if( ( data_number = ( int * * ) malloc( data_row * sizeof( int * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < data_row; index++ )
    {
	if( ( *( data_number + index ) = ( int * ) malloc( 2 * sizeof( int ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    if( ( data_src = ( double * * ) malloc( data_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < data_row; index++ )
    {
	if( ( *( data_src + index ) = ( double * ) malloc( data_column * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // graph traversal
    for( int index = 0, index_temp = 0; index < graph->count_vertex && index_temp < data_row; index++ )
    {
	AdjListNode * head_node = graph->array[ index ].head;

	while( head_node != NULL )
	{
	    // data_number
	    data_number[ index_temp ][ 0 ] = index;
	    data_number[ index_temp ][ 1 ] = head_node->dest;

	    // data_src
	    for( int index_i = 0; index_i < data_column; index_i++ )
	    {
		data_src[ index_temp ][ index_i ] = head_node->weight[ index_i ];
	    }

	    head_node = head_node->next;
	    index_temp++;
	}
    }
#if 1    // print data
    puts( ">>>>>>>>>>>>data_number:" );
    for( int index_i = 0; index_i < data_row; index_i++ )
    {
	for( int index_j = 0; index_j < 2; index_j++ )
	{
	    printf( "%d ", data_number[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
    puts( "\n>>>>>>>>>>>>data_src:" );
    for( int index_i = 0; index_i < data_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_column; index_j++ )
	{
	    printf( "%.4le ", data_src[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif

    adjustment_impl( data_number, data_src, number_known, data_known,
	    data_row, data_column, count_known, graph->count_vertex);

    // free
    for( int index = 0; index < data_row; index++ )
    {
	free( *( data_src + index ) );
    }
    free( data_src );
    for( int index = 0; index < data_row; index++ )
    {
	free( *( data_number + index ) );
    }
    free( data_number );
}

void adjustment_impl( int * * data_number, double * * data_src, int * number_known, double * * data_known,
	int data_row, int data_column, int count_known, int count )
{
    puts( "============adjustment of GNSS network implementation============" );

    /*
     * count of baseline observation data contains unknown node(s)
     * */
    int row_count_unknown = 0;
    for( int index_i = 0; index_i < data_row; index_i++ )
    {
	//if( ( int status = IsKnownVertex( *( data_number + index_i ), number_known, 2, count_known ) ) != 1 )
	if( IsKnownVertex( *( data_number + index_i ), number_known, 2, count_known ) != 1 )
	{
	    row_count_unknown++;
	}
    }
#if 1
    printf( ">>>>row count of unknown is %d\n", row_count_unknown );
#endif

    // linear system
    /*
     * count_node_unknwon = count - count_known
     *
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
    // int data_mat_row = row_count_unknown * 3, data_mat_column = ( count - count_known ) * 3;
    int data_mat_row = data_row * 3, data_mat_column = ( count - count_known ) * 3;
    double * * data_mat = NULL;
    double * * data_weight = NULL;
    double * data_rhs = NULL;
    double * data_sol = NULL;

    // store original linear system
    /*
     * residual = original_rhs - original_mat * data_sol
     * */
    double * * original_mat = NULL;    // size data_mat_row x data_mat_column
    double * original_rhs = NULL;    // size data_mat_row x 1
    double * original_res = NULL;    // size data_mat_row x 1

    // memory allocation for original_rhs, original_res
    if( ( original_rhs = ( double * ) malloc( data_mat_row * sizeof( double ) ) ) 
	    == NULL 
	    || ( original_res = ( double * ) malloc( data_mat_row * sizeof( double ) ) ) == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize original_rhs = 0, original_res = 0
    for( int index = 0; index < data_mat_row; index++ )
    {
	original_rhs[ index ] = 0;
	original_res[ index ] = 0;
    }

    if( ( original_mat = ( double * * ) malloc( data_mat_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < data_mat_row; index++ )
    {
	if( ( *( original_mat + index ) = ( double * ) malloc( data_mat_column * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocatio failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // initialize original_mat = 0
    for( int index_i = 0; index_i < data_mat_row; index_i++ )
    {
	for( int index_j = 0; index_j < data_mat_column; index_j++ )
	{
	    original_mat[ index_i ][ index_j ] = 0;
	}
    }

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
#if 0    // print data_mat
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
#if 0    // print data_weight
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
#if 0    // print data_rhs
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
#if 1
    assemble_adjust_linsys( data_mat, data_rhs, data_weight,
	    data_number, data_src,
	    number_known, data_known,
	    data_row, count_known, count, row_count_unknown,
	    original_mat, original_rhs );
#endif

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
#if 1
    solver_adjust_linsys( data_mat, data_rhs, data_mat_row, data_mat_column,
	    data_sol );
#endif

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

    // computing residual
    /*
     * original_res = original_rhs - original_mat * data_sol
     * */
    computing_residual( original_res, original_rhs, original_mat, data_sol,
	    data_mat_row, data_mat_column );
#if 1    // print residual
    puts( "\n>>>>>>>>>>>>residual:" );
    for( int index = 0; index < data_mat_row; index++ )
    {
	printf( "%.4le\n", original_res[ index ] );
    }
#endif

    // free
    for( int index = 0; index < data_row; index++ )
    {
	free( *( original_mat + index ) );
    }
    free( original_mat );
    free( original_rhs );
    free( original_res );
    for( int index = 0; index <data_mat_row; index++ )
    {
	free( *( data_weight + index ) );
    }
    free( data_weight );
    free( data_sol );
    free( data_rhs );
    for( int index = 0; index < data_mat_row; index++ )
    {
	free( *( data_mat + index ) );
    }
    free( data_mat );
}

int IsKnownVertex( int * vertex, int * known_vertex, int size_vertex, int size_known_vertex )
{
    int status = 1;
    int index_i = 0, index_j = 0;

    for( index_i = 0; index_i < size_vertex; index_i++ )
    {
	for( index_j = 0; index_j < size_known_vertex; index_j++ )
	{
	    if( *( known_vertex + index_j ) == *( vertex + index_i ) )
	    {
		break;
	    }
	}

	if( index_j == size_known_vertex )
	{
	    status = 0;
	}

    }

    return status;
}
