#include "../include/main.h"

int main()
{
    puts( "============test adjustment of GNSS network==========" );

    int count = 0;    // count of vertices
    printf( ">>>>>>>>>>>>input count of vertices: " );
    scanf( "%d", &count );
    AdjGraph * value_graph = GraphGeneration( count );    // graph generation with count vertices

    FILE * fp_adjust = NULL;    // observation data file
    if( ( fp_adjust = fopen( "../file/adjust_data.txt", "rb"  ) ) 
    //if( ( fp_adjust = fopen( "../file/adjust_data_1.txt", "rb"  ) ) 
	    == NULL )
    {
	fprintf( stderr, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    int size_data_row = 0, size_data_column = 0;
    fscanf( fp_adjust, "%d%d", &size_data_row, &size_data_column );

#if 1
    printf( "row = %d\ncolumn = %d\n", size_data_row, size_data_column );
#endif

    int node_location[ 2 ] = { 0 };    // node enumeration
    double weight_data[ 9 ] = { 0 };    // weight data in edge
    for( int index = 0; index < size_data_row; index++ )
    {
	for( int index_i = 0; index_i < 2; index_i++ )
	{
	    fscanf( fp_adjust, "%d", node_location + index_i );
	}
	for( int index_i = 0; index_i < size_data_column - 2; index_i++ )
	{
	    fscanf( fp_adjust, "%lf", weight_data + index_i );
	}

	( value_graph->count_edge )++;    // count of edge augment by 1
        GraphInsert( value_graph, node_location, weight_data, 9 );
    }

    // fclose
    fclose( fp_adjust );

#if 0    // graph display
    puts( ">>>>>>>>>>>>graph display:" );
    GraphDisplay( value_graph );
#endif

    int size_known_vertex = 0;    // count of known vertex
    printf( "\n>>>>>>>>>>>>input count of known vertex: " );
    scanf( "%d", &size_known_vertex );

    int * number_known_vertex = NULL;    // number of known vertex
    double * * coordinate_known_vertex = NULL;    // coordinate of known vertex

    /*
     * size number_known_vertex = size_known_vertex x 1
     * */
    if( ( number_known_vertex = ( int * ) malloc( size_known_vertex * sizeof( int ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!( number_known_vertex )\n" );
	exit( EXIT_FAILURE );
    }

    /*
     * size coordinate_known_vertex = size_known_vertex x 3
     * */
    if( ( coordinate_known_vertex = ( double * * ) malloc( size_known_vertex * sizeof( double * ) ) ) == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < size_known_vertex; index++ )
    {
	if( ( *( coordinate_known_vertex + index ) = ( double * ) malloc( 3 * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    printf( "\n>>>>>>>>>>>>input information of known vertex:\n" );
    for( int index = 0; index < size_known_vertex; index++ )
    {
	printf( "number of known vertex: " );
	scanf( "%d", number_known_vertex + index );

	printf( "\ncoordinate of known vertex: " );
	for( int index_i = 0; index_i < 3; index_i++ )
	{
	    scanf( "%lf", *( coordinate_known_vertex + index ) + index_i );
	}
    }

    graph_adjustment( value_graph, size_known_vertex, number_known_vertex, coordinate_known_vertex );

    // free, fclose
    for( int index = 0; index < size_known_vertex; index++ )
    {
	free( *( coordinate_known_vertex + index ) );
    }
    free( coordinate_known_vertex );
    free( number_known_vertex );

    return 0;
}
