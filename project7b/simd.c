/***********************************************************
 * Program: simd.c
 * Author: Joel Huffman
 * Last updated: 6/6/2019
 * Sources: http://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj07b.html
 ***********************************************************/
#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <iostream>
#include "simd.p4.h"

#define NUMTRIES 20

#ifndef _OPENMP
	fprintf( stderr, "No OpenMP support!\n" );
	exit (1);
#endif

int main()
{
	int     Size;
	float * Array;
	float * Sums;
	FILE *  fp;
	int     i;

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

	// run SIMD work

	double maxPerformance = 0.;
	double avgPerformance = 0.;

	// run multiple times to get average performance as well as peak performance
	for (int t = 0; t < NUMTRIES; t++)
	{
		// start time
		double time0 = omp_get_wtime();

		for (int shift = 0; shift < Size; shift++)
		{
			Sums[shift] = SimdMulSum(Array, &Array[shift], Size);
		}

		// end time
		double time1 = omp_get_wtime();

		double megaCalcsPerSecond = (double) (Size * Size) / (time1 - time0) / 1000000.;
		if (megaCalcsPerSecond > maxPerformance)
		{
			maxPerformance = megaCalcsPerSecond;
		}
		avgPerformance += megaCalcsPerSecond;
	}

	avgPerformance /= NUMTRIES;

	// print results
	printf ("%lf\t%lf\n",avgPerformance, maxPerformance);

	// append results to results.txt
	FILE* resultsFp = fopen("results.txt", "a");
	fprintf (resultsFp, "%lf\t%lf\n", avgPerformance, maxPerformance);
	fclose(resultsFp);
}
