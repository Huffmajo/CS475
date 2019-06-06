/***********************************************************
 * Program: openmp.c
 * Author: Joel Huffman
 * Last updated: 6/6/2019
 * Sources: http://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj07b.html
 ***********************************************************/
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>

#define NUMTRIES 20

#ifndef _OPENMP
	fprintf( stderr, "No OpenMP support!\n" );
	return 1;
#endif

int main ()
{
	int     Size;
	float * Array;
	float * Sums;
	FILE *  fp;
	int     i;

	// set number of threads to utilize
	omp_set_num_threads(NUMT);

	// read in and store values from signal.txt
	fp = fopen( "signal.txt", "r" );
	if( fp == NULL )
	{
		fprintf( stderr, "Cannot open file 'signal.txt'\n" );
		exit( 1 );
	}
	fscanf( fp, "%d", &Size );
	Array = (float *)malloc( 2 * Size * sizeof(float) );
	Sums  = (float *)malloc( 1 * Size * sizeof(float) );
	for( i = 0; i < Size; i++ )
	{
		fscanf( fp, "%f", &Array[i] );
		Array[i+Size] = Array[i];		// duplicate the array
	}
	fclose( fp );

	// do openMP work

	double maxPerformance = 0.;
	double avgPerformance = 0.;

	// run multiple times to get average performance as well as peak performance
	for (int t = 0; t < NUMTRIES; t++)
	{
		// start time
		double time0 = omp_get_wtime();

		#pragma omp parallel for default(none) shared(Size, Array, Sums)
		for (int shift = 0; shift < Size; shift++)
		{
			float sum = 0.;
			for (int i = 0; i < Size; i++)
			{
				sum += Array[i] * Array[i + shift];
			}

			Sums[shift] = sum;
		}

		// end time
		double time1 = omp_get_wtime();

		double megaCalcsPerSecond = (double)Size / (time1 - time0) / 1000000.;
		if (megaCalcsPerSecond > maxPerformance)
		{
			maxPerformance = megaCalcsPerSecond;
		}
		avgPerformance += megaCalcsPerSecond;
	}

	// print Sums[x] vs shift data to monitor for sine wave
	if (NUMT == 32)
	{
		printf("Autocorrelate Data\n");
		printf("Sums[x] Shift\n");

		FILE* autoFp = fopen("autocorrelateResults.txt", "a");
		for (int i = 1; i < 513; i++)
		{
			printf ("%d\t%lf\n", i, Sums[i]);
			fprintf (autoFp, "%d\t%lf\n", i, Sums[i]);
		}
		fclose(autoFp);
	}

	avgPerformance /= NUMTRIES;

	// print results
	printf ("%d\t%lf\t%lf\n", NUMT, avgPerformance, maxPerformance);

	// append results to results.txt
	FILE* resultsFp = fopen("results.txt", "a");
	fprintf (resultsFp, "%d\t%lf\t%lf\n", NUMT, avgPerformance, maxPerformance);
	fclose(resultsFp);

	return 0;
}
