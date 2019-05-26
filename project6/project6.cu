/***********************************************************
 * Program: project6.cu
 * Author: Joel Huffman
 * Last updated: 5/25/2019
 * Sources: http://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj06.html
 ***********************************************************/

// System includes
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include "helper_functions.h"
#include "helper_cuda.h"


#ifndef BLOCKSIZE
#define BLOCKSIZE		32		// number of threads per block
#endif

#ifndef SIZE
#define SIZE			1*1024*1024	// array size
#endif

#ifndef NUMTRIALS
#define NUMTRIALS		100		// to make the timing more accurate
#endif

#ifndef TOLERANCE
#define TOLERANCE		0.00001f	// tolerance to relative error
#endif

// ranges for the random numbers:
const float XCMIN =	 0.0;
const float XCMAX =	 2.0;
const float YCMIN =	 0.0;
const float YCMAX =	 2.0;
const float RMIN  =	 0.5;
const float RMAX  =	 2.0;

float Ranf( float low, float high )
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}

int Ranf( int ilow, int ihigh )
{
        float low = (float)ilow;
        float high = ceil( (float)ihigh );

        return (int) Ranf(low,high);
}

void TimeOfDaySeed( )
{
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  timer;
	time( &timer );
	double seconds = difftime( timer, mktime(&y2k) );
	unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
	srand( seed );
}

// determine if vertical plate is hit by laser
__global__  void MonteCarlo( float *A, float *B, float *C, float *D )
{
/*
	__shared__ float numHits[BLOCKSIZE];
	unsigned int numItems = blockDim.x;
	unsigned int tnum = threadIdx.x;
	unsigned int wgNum = blockIdx.x;
*/
	unsigned int gid = blockIdx.x*blockDim.x + threadIdx.x;

	// randomize the location and radius of the circle:
	float xc = A[gid];
	float yc = B[gid];
	float r =  C[gid];

	// solve for the intersection using the quadratic formula:
	float a = 2.;
	float b = -2.*( xc + yc );
	float c = xc*xc + yc*yc - r*r;
	float d = b*b - 4.*a*c;

	// If d is less than 0, then the circle was completely missed (Case A) 
	if (d >= 0.)
	{
		// hits the circle:
		// get the first intersection:
		d = sqrtf( d );
		float t1 = (-b + d ) / ( 2.*a );	// time to intersect the circle
		float t2 = (-b - d ) / ( 2.*a );	// time to intersect the circle
		float tmin = t1 < t2 ? t1 : t2;		// only care about the first intersection

		// If tmin is less than 0., then the circle completely engulfs the laser pointer (Case B)
		if (tmin >= 0.)
		{
			// where does it intersect the circle?
			float xcir = tmin;
			float ycir = tmin;

			// get the unitized normal vector at the point of intersection:
			float nx = xcir - xc;
			float ny = ycir - yc;
			float n = sqrtf( nx*nx + ny*ny );
			nx /= n;	// unit vector
			ny /= n;	// unit vector

			// get the unitized incoming vector:
			float inx = xcir - 0.;
			float iny = ycir - 0.;
			float in = sqrtf( inx*inx + iny*iny );
			inx /= in;	// unit vector
			iny /= in;	// unit vector

			// get the outgoing (bounced) vector:
			float dot = inx*nx + iny*ny;
//			float outx = inx - 2.*nx*dot;	// angle of reflection = angle of incidence`
			float outy = iny - 2.*ny*dot;	// angle of reflection = angle of incidence`

			// find out if it hits the infinite plate:
			float t = ( 0. - ycir ) / outy;

			// If t is less than 0., then the reflected beam went up instead of down (Case C)
			if (t >= 0.)
			{
				D[gid] = 1;;
			}
		}
	}

/*
	prods[tnum] = A[gid] * B[gid];

	for (int offset = 1; offset < numItems; offset *= 2)
	{
		int mask = 2 * offset - 1;
		__syncthreads();
		if ((tnum & mask) == 0)
		{
			numHits[tnum] += numHits[tnum + offset];
		}
	}

	__syncthreads();
	if (tnum == 0)
		D[wgNum] = numHits[0];
*/
}


// main program:

int
main( int argc, char* argv[ ] )
{
//	int dev = findCudaDevice(argc, (const char **)argv);

	// allocate host memory:
	float *xcs = new float [ NUMTRIALS ];
	float *ycs = new float [ NUMTRIALS ];
	float *rs = new float [ NUMTRIALS ];
	float *hits = new float [ NUMTRIALS ];

	// fill arrays with random values in range
	for( int n = 0; n < NUMTRIALS; n++ )
	{
		xcs[n] = Ranf( XCMIN, XCMAX );
                ycs[n] = Ranf( YCMIN, YCMAX );
                rs[n] = Ranf(  RMIN,  RMAX ); 
		hits[n] = 0.;
	}

	// allocate device memory:

	float *dxcs, *dycs, *drs, *dhits;

	dim3 dimsxcs( NUMTRIALS, 1, 1 );
	dim3 dimsycs( NUMTRIALS, 1, 1 );
	dim3 dimsrc( NUMTRIALS, 1, 1 );
	dim3 dimshits( NUMTRIALS, 1, 1 );

	//__shared__ float prods[SIZE/BLOCKSIZE];


	cudaError_t status;
	status = cudaMalloc( reinterpret_cast<void **>(&dxcs), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&dycs), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&drs), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&dhits), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );

	// copy host memory to the device:

	status = cudaMemcpy( dxcs, xcs, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( dycs, ycs, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( drs, rs, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( dhits, hits, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );

	// setup the execution parameters:

	dim3 threads(BLOCKSIZE, 1, 1 );
	dim3 grid( NUMTRIALS / threads.x, 1, 1 );

	// Create and start timer

	cudaDeviceSynchronize( );

	// allocate CUDA events that we'll use for timing:

	cudaEvent_t start, stop;
	status = cudaEventCreate( &start );
		checkCudaErrors( status );
	status = cudaEventCreate( &stop );
		checkCudaErrors( status );

	// record the start event:

	status = cudaEventRecord( start, NULL );
		checkCudaErrors( status );

	// execute the kernel:

//	for( int t = 0; t < NUMTRIALS; t++)
//	{
	        MonteCarlo<<< grid, threads >>>( dxcs, dycs, drs, dhits );
//	}

	// record the stop event:

	status = cudaEventRecord( stop, NULL );
		checkCudaErrors( status );

	// wait for the stop event to complete:

	status = cudaEventSynchronize( stop );
		checkCudaErrors( status );

	float msecTotal = 0.0f;
	status = cudaEventElapsedTime( &msecTotal, start, stop );
		checkCudaErrors( status );

	// copy result from the device to the host:

	status = cudaMemcpy( hits, dhits, NUMTRIALS*sizeof(float), cudaMemcpyDeviceToHost );
		checkCudaErrors( status );

	// add up all the hits
	int numHits = 0;
	for (int i = 0; i < NUMTRIALS; i++)
	{
		if (hits[i] == 1)
		{
			numHits++;
		}
	}

	// compute and print the performance
	double secondsTotal = 0.001 * (double)msecTotal;
	double multsPerSecond = (float)NUMTRIALS / secondsTotal;
	double megaMultsPerSecond = multsPerSecond / 1000000.;
	double probability = (float)numHits / (float)NUMTRIALS;

	// print performance
	printf("%d\t%d\t%lf\t%lf\n", NUMTRIALS, BLOCKSIZE, megaMultsPerSecond, probability);

	// also write performance to results.txt
	FILE *fp;
	fp = fopen("results.txt", "a");
	fprintf(fp, "%d\t%d\t%lf\t%lf\n", NUMTRIALS, BLOCKSIZE, megaMultsPerSecond, probability);
	fclose(fp);

	// clean up memory:
	delete [ ] xcs;
	delete [ ] ycs;
	delete [ ] rs;
	delete [ ] hits;

	status = cudaFree( dxcs );
		checkCudaErrors( status );
	status = cudaFree( dycs );
		checkCudaErrors( status );
	status = cudaFree( drs );
		checkCudaErrors( status );
	status = cudaFree( dhits );
		checkCudaErrors( status );

	return 0;
}

