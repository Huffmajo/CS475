/***********************************************************
 * Program: project0.c
 * Author: Joel Huffman
 * Last updated: 4/1/2019
 * Sources: https://oregonstate.instructure.com/courses/1716463/pages/week-01-listen-to-this?module_item_id=18618494
 * https://www.gamedev.net/forums/topic/392211-max-value-for-double-/
 ***********************************************************/
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define NUMTRIES 200
#define ARRAYSIZE 10000000

float A[ARRAYSIZE];
float B[ARRAYSIZE];
float C[ARRAYSIZE];

int
main()
{
#ifndef _OPENMP
	fprintf(stderr, "OpenMP is not supported here -- sorry.\n");
	return 1;
#endif

	// find out how many cores are on system
	printf("\nNumber of system cores: %d\n", omp_get_num_procs());

	omp_set_num_threads(NUMT);
	printf("Using %d threads with sample size of %d\n", NUMT, NUMTRIES);

	double maxMegaMults = 0.;
	double sumMegaMults = 0.;
	double fastestTime = DBL_MAX;
	double sumFastestTime = 0.;

	for (int t = 0; t < NUMTRIES; t++)
	{
		double time0 = omp_get_wtime();

#pragma omp parallel for
		for (int i = 0; i < ARRAYSIZE; i++)
		{
			C[i] = A[i] * B[i];
		}

		double time1 = omp_get_wtime();
		double megaMults = (double)ARRAYSIZE / (time1 - time0) / 1000000.;
		double elapsedTime = (time1 - time0) * 1000000;
		
		if (megaMults > maxMegaMults)
		{
			maxMegaMults = megaMults;
		}

		if (elapsedTime < fastestTime)
		{
			fastestTime = elapsedTime;
		}

		sumMegaMults += megaMults;
		sumFastestTime += fastestTime;

	}

	printf("Peak Performance = %8.2lf MegaMults/Sec\n", maxMegaMults);
	printf("Average Performance = %8.2lf MegaMults/Sec\n", (sumMegaMults / (double) NUMTRIES));
	printf("Fastest Time = %8.2lf μsec\n", fastestTime);
	printf("Average Time = %8.2lf μsec\n", (sumFastestTime / (double) NUMTRIES));
	return 0;
}
