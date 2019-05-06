/***********************************************************
 * Program: project3.c
 * Author: Joel Huffman
 * Last updated: 5/05/2019
 * Sources: http://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj03.html
 ***********************************************************/
#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>

#define NUMT 4

// global variables
int	NowYear;		// 2019 - 2024
int	NowMonth;		// 0 - 11
float	NowPrecip;		// inches of rain per month
float	NowTemp;		// temperature this month
float	NowHeight;		// grain height in inches
int	NowNumDeer;		// number of deer in the current population
int	NowNumSasquatch;	// number of sasquatch in the current population

const float GRAIN_GROWS_PER_MONTH =		8.0;
const float ONE_DEER_EATS_PER_MONTH =		0.5;

const float AVG_PRECIP_PER_MONTH =		6.0;	// average
const float AMP_PRECIP_PER_MONTH =		6.0;	// plus or minus
const float RANDOM_PRECIP =			2.0;	// plus or minus noise

const float AVG_TEMP =				50.0;	// average
const float AMP_TEMP =				20.0;	// plus or minus
const float RANDOM_TEMP =			10.0;	// plus or minus noise

const float MIDTEMP =				40.0;
const float MIDPRECIP =				10.0;

const float ONE_SASQUATCH_EATS_PER_MONTH =	1.5;

// returns squared result of provided float
float SQR( float x )
{
        return x*x;
}

// returns random float within low to high range
float RandFloat( unsigned int *seedp,  float low, float high )
{
        float r = (float) rand_r( seedp );              // 0 - RAND_MAX

        return(low + r * (high - low) / (float)RAND_MAX);
}

// returns random integer within low to high range
/*
int RandInt( unsigned int *seedp, int ilow, int ihigh )
{
        float low = (float)ilow;
        float high = (float)ihigh + 0.9999f;

        return (int)(RandFloat(seedp, low,high) );
}
*/

// manages deer thread
void GrainDeer()
{
	while( NowYear < 2025 )
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:

		int nextNumDeer;
		float nextHeight;

		// if more graindeer than inches of grain, graindeer population decrements
		if (NowNumDeer > NowHeight)
		{
			nextNumDeer = NowNumDeer - 1;
		}

		// if fewer graindeer than inches of grain, graindeer population increments
		else if (NowNumDeer < NowHeight)
		{
			nextNumDeer = NowNumDeer + 1;
			
			// we can't have negative deer
			if (nextNumDeer < 0)
			{
				nextNumDeer = 0;
			}
		}

		// DoneComputing barrier:
		#pragma omp barrier

		// copy next state into now state
		NowNumDeer = nextNumDeer;

		// DoneAssigning barrier:
		#pragma omp barrier

		// nothing needed here

		// DonePrinting barrier:
		#pragma omp barrier

		// nothing needed here
	}
}

// manages grain thread
void Grain()
{
	while( NowYear < 2025 )
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:

		float nextHeight;
		float tempFactor = exp(   -SQR(  ( NowTemp - MIDTEMP ) / 10.  )   );
		float precipFactor = exp(   -SQR(  ( NowPrecip - MIDPRECIP ) / 10.  )   );

		// compute next grain height
		nextHeight = NowHeight + (tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH);
		nextHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
		nextHeight -= (float)NowNumSasquatch * ONE_SASQUATCH_EATS_PER_MONTH;

		// can't have negative height
		if (nextHeight < 0.)
		{
			nextHeight = 0.;
		}

		// DoneComputing barrier:
		#pragma omp barrier

		// copy next state into now state
		NowHeight = nextHeight;

		// DoneAssigning barrier:
		#pragma omp barrier

		// nothing needed here

		// DonePrinting barrier:
		#pragma omp barrier

		// nothing needed here
	}
}

// updates timing and prints results of other threads
void Watcher()
{
	unsigned int seed = 0;
	while( NowYear < 2025 )
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:

		// nothing needed here		

		// DoneComputing barrier:
		#pragma omp barrier

		// nothing needed here

		// DoneAssigning barrier:
		#pragma omp barrier

		// calculate and update weather
		float ang = (  30.*(float)NowMonth + 15.  ) * ( M_PI / 180. );

		float temp = AVG_TEMP - AMP_TEMP * cos( ang );
		NowTemp = temp + RandFloat( &seed, -RANDOM_TEMP, RANDOM_TEMP );

		float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang );
		NowPrecip = precip + RandFloat( &seed,  -RANDOM_PRECIP, RANDOM_PRECIP );
		if( NowPrecip < 0. )
			NowPrecip = 0.;	

		//convert to metric for printing
		float metricHeight = NowHeight * 2.54;
		float metricPrecip = NowPrecip * 2.54;
		float metricTemp = (NowTemp - 32) * (5. / 9.);

		// print results for this month
		printf("%d/%d\t%d\t%d\t%lf\t%lf\t%lf\n", NowMonth + 1, NowYear, NowNumDeer, NowNumSasquatch, metricHeight, metricPrecip, metricTemp);

		// write results to .txt file as well
		FILE *fp;
		fp = fopen("results.txt", "a");
		fprintf(fp, "%d/%d\t%d\t%d\t%lf\t%lf\t%lf\n", NowMonth + 1, NowYear, NowNumDeer, NowNumSasquatch, metricHeight, metricPrecip, metricTemp);
		fclose(fp);

		// increment month
		NowMonth++;

		// account for last month of the year
		if (NowMonth > 11)
		{
			NowYear++;
			NowMonth = 0;
		}

		// DonePrinting barrier:
		#pragma omp barrier

		// nothing needed here
	}
}

// manages sasquatch (custom agent) thread
void Sasquatch()
{
	while( NowYear < 2025 )
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:

		int nextNumSasquatch;
		
		// if the deer-to-sasquatch ratio is over 2:1, increment sasquatch
		if (NowNumDeer > NowNumSasquatch * 2)
		{
			nextNumSasquatch = NowNumSasquatch + 1;
		}

		// if there are fewer than 2:1 deer-to-sasquatch, remove all sasquatch
		else if (NowNumDeer < NowNumSasquatch * 2)
		{
			nextNumSasquatch = 0;

			if (nextNumSasquatch < 0)
			{
				nextNumSasquatch = 0;
			}
		}

		// DoneComputing barrier:
		#pragma omp barrier

		// copy next state into now state
		NowNumSasquatch = nextNumSasquatch;

		// DoneAssigning barrier:
		#pragma omp barrier

		// nothing needed here

		// DonePrinting barrier:
		#pragma omp barrier

		// nothing needed here
	}
}

int main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
	fprintf( stderr, "No OpenMP support!\n" );
	return 1;
#endif

	//set number of threads
	omp_set_num_threads(NUMT);

	// starting date and time:
	NowMonth =    0;
	NowYear  = 2019;

	// starting state (feel free to change this if you want):
	NowNumDeer = 2;
	NowHeight =  1.5;


	#pragma omp parallel sections
	{
		#pragma omp section
		{
			GrainDeer( );
		}

		#pragma omp section
		{
			Grain( );
		}

		#pragma omp section
		{
			Watcher( );
		}

		#pragma omp section
		{
			Sasquatch( );	
		}

	}       // implied barrier -- all functions must return in order
		// to allow any of them to get past here
}
