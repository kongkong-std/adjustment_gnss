#include <stdio.h>
#include <stdlib.h>

/*
 * observation data
 * count of observation data
 * every row contains x, y, z, hessian matrix
 * known point
 * */
void adjustment_impl( int **, double **, int, int, int );

int main()
{
    puts( "============test adjustment of GNSS network==========" );

    FILE * fp_adjust = NULL;
    if( ( fp_adjust = fopen( "../file/adjust_data.txt", "rb"  ) ) 
	    == NULL )
    {
	fprintf( stderr, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    int size_data_row = 0, size_data_column = 0;
    fscanf( fp_adjust, "%d%d", &size_data_row, &size_data_column );

#if 0
    printf( "row = %d\ncolumn = %d\n", size_data_row, size_data_column );
#endif

    int * * data_number = NULL;    // observation number
    double * * data_adjust = NULL;    // observation data

    // data_number allocation
    /*
     * size data_number = size_data_row x 2
     * */
    if( ( data_number = ( int * * ) malloc( size_data_row * sizeof( int * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < size_data_row; index++ )
    {
	if( ( *( data_number + index ) = ( int * ) malloc( 2 * sizeof( int ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // data_adjust memory allocation
    /*
     * size data_adjust = size_data_row x ( size_data_column - 2 )
     * */
    if( ( data_adjust = ( double * * ) malloc( size_data_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < size_data_row; index++ )
    {
	if( ( *( data_adjust + index ) = ( double * ) malloc( ( size_data_column - 2 ) * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    for( int index = 0; index < size_data_row; index++ )
    {
	for( int index_i = 0; index_i < 2; index_i++ )
	{
	    fscanf( fp_adjust, "%d", *( data_number + index ) + index_i );
	}
	for( int index_j = 0; index_j < size_data_column - 2; index_j++ )
	{
	    fscanf( fp_adjust, "%lf", *( data_adjust + index ) + index_j );
	}
    }

#if 0    // print data_src
    puts( ">>>>show source data" );
    for( int index = 0; index < size_data_row; index++ )
    {
	for( int index_i = 0; index_i < 2; index_i++ )
	{
	    printf( "%d ", data_number[ index ][ index_i ] );
	}
	for( int index_j = 0; index_j < size_data_column - 2; index_j++ )
	{
	    printf( "%lf ", data_adjust[ index ][ index_j ] );
	}
	putchar( '\n' );
    }
#endif

    int size_data_known = 2;    // count of known point

    adjustment_impl( data_number, data_adjust, size_data_row, size_data_column, size_data_known );

    // free, fclose
    for( int index = 0; index < size_data_row; index++ )
    {
	free( *( data_number + index ) );
    }
    free( data_number );
    for( int index = 0; index < size_data_row; index++ )
    {
	free( *( data_adjust + index ) );
    }
    free( data_adjust );
    fclose( fp_adjust );

    return 0;
}
