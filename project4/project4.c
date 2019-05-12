/***********************************************************
 * Program: project4.c
 * Author: Joel Huffman
 * Last updated: 5/11/2019
 * Sources: http://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj04.html
 ***********************************************************/
#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include "simd.p4.h"

#define NUMTRIES 50

void MultNonSIMD(int size, float* a, float* b, float* c)
{
	float sum = 0.;

	for (int i = 0; i < size; i++)
	{
		// array multiplication
		c[i] = a[i] * b[i];
	}	
}

float MultReductNonSIMD(int size, float* a, float* b)
{
	float sum = 0.;

	for (int i = 0; i < size; i++)
	{
		// array multiplication
		sum += a[i] * b[i];
	}

	return sum;
}



int main( int argc, char *argv[ ] )
{
	#ifndef _OPENMP
		fprintf( stderr, "No OpenMP support!\n" );
		return 1;
	#endif

	double peakMegaMults = 0.;
	double avgMegaMults = 0.;
	float* a = new float [ARRAY_SIZE];
	float* b = new float [ARRAY_SIZE];
	float* c = new float [ARRAY_SIZE];

	// open results.txt to append results to later
	FILE *fp;
	fp = fopen("results.txt", "a");

	// run multiple times to determine average and peak
	for (int i = 0; i < NUMTRIES; i++)
	{
		// start timer
		double startTime = omp_get_wtime();

		// run the correct function for the method we are testing
		switch(METHOD)
		{
			// SIMD multiply
			case 0:
				SimdMul(a, b, c, ARRAY_SIZE);
				break;

			// non-SIMD multiply
			case 1:
				MultNonSIMD(ARRAY_SIZE, a, b, c);
				break;

			// SIMD multiply with reduction
			case 2:
				SimdMulSum(a, b, ARRAY_SIZE);
				break;

			// non-SIMD multiply with reduction
			case 3:
				MultReductNonSIMD(ARRAY_SIZE, a, b);
				break;

			// if METHOD is an unexpected number
			default:
				printf("Something has gone wrong. Unknown method number %d called\n", METHOD);
				exit(1);
		}

		// stop timer
		double endTime = omp_get_wtime();

		// get time and performance of run
		double megaMults = (double)ARRAY_SIZE / (endTime - startTime) / 1000000.;

		avgMegaMults += megaMults;

		if (megaMults > peakMegaMults)
		{
			peakMegaMults = megaMults;
		}
	}
	
	// calculate average to gauge validity
	avgMegaMults /= (double)NUMTRIES;

	// print results
	printf("%d\t%lf\t%lf\n", ARRAY_SIZE, avgMegaMults, peakMegaMults);

	// append results to results.txt
	fprintf(fp, "%d\t%lf\t%lf\n", ARRAY_SIZE, avgMegaMults, peakMegaMults);
	fclose(fp);
}
