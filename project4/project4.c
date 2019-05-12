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

void MultNonSIMD(int size, double* a, double* b, double* c, int reductFlag)
{
	// if reductFlag is 1, apply reduction to method
	

}



int main( int argc, char *argv[ ] )
{
	#ifndef _OPENMP
		fprintf( stderr, "No OpenMP support!\n" );
		return 1;
	#endif

	double peakMegaMults = 0.;
	double avgMegaMults = 0.;
	double* a[ARRAY_SIZE];
	double* b[ARRAY_SIZE];
	double* c[ARRAY_SIZE];

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
				MultNonSIMD(ARRAY_SIZE, a, b, c, 0);
				break;

			// SIMD multiply with reduction
			case 2:
				SimdMulSum(a, b, c, ARRAY_SIZE);
				break;

			// non-SIMD multiply with reduction
			case 3:
				MultNonSIMD(ARRAY_SIZE, a, b, c, 1);
				break;

			// if METHOD is an unexpected number
			default:
				printf("Something has gone wrong. Unknown method number %d called\n", METHOD);
				exit(1);
				break;
		}

/*
		if (METHOD == 0)
		{

		}
		
		else if ()
		{

		}

		else if ()
		{

		}

		else if ()
		{

		}

		else 
		{
			printf("Unknown method %d called\n", METHOD);
			exit 1;
		}
*/
		// stop timer
		double endTime = omp_get_wtime();
	}
}
