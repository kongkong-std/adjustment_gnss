#include <stdio.h>
#include <stdlib.h>

void adjustment_implementation( double **, int, int );

int main()
{
    puts( "========test adjustment of GNSS network========" );

    FILE * fp_adjust = NULL;
    if( ( fp_adjust = fopen( "../file/adjustment_data.txt", "rb" ) ) 
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

    double ** adjust_data = NULL;
    if( ( adjust_data = ( double * * ) malloc( size_data_row * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( int index = 0; index < size_data_row; index++ )
    {
	if( ( *( adjust_data + index ) = ( double * ) malloc( size_data_column * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    for( int index_i = 0; index_i < size_data_row; index_i++ )
    {
	for( int index_j = 0; index_j < size_data_column; index_j++ )
	{
	    fscanf( fp_adjust, "%lf", *( adjust_data + index_i ) + index_j );
	}
    }

#if 0
    for( int index_i = 0; index_i < size_data_row; index_i++ )
    {
	for( int index_j = 0; index_j < size_data_column; index_j++ )
	{
	    printf( "%.4lf ", *( *( adjust_data + index_i ) + index_j ) );
	}
	putchar( '\n' );
    }
#endif

    adjustment_implementation( adjust_data, size_data_row, size_data_column );

    // free and fclose
    for( int index = 0; index < size_data_row; index++ )
    {
	free( *( adjust_data + index ) );
    }
    free( adjust_data );
    fclose( fp_adjust );

    return 0;
}
