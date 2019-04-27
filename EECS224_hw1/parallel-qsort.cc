/**
 *  \file parallel-qsort.cc
 *
 *  \brief Implement your parallel quicksort algorithm in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include <omp.h>

#include "sort.hh"


#define NUM_PROC 8

/**
 *   Given a pivot value, this routine partitions a given input array
 *   into two sets: the set A_le, which consists of all elements less
 *   than or equal to the pivot, and the set A_gt, which consists of
 *   all elements strictly greater than the pivot.
 *      
 *   This routine overwrites the original input array with the
 *   partitioned output. It also returns the index n_le such that
 *   (A[0:(k-1)] == A_le) and (A[k:(N-1)] == A_gt).
 */

int paraPartition(keytype pivot, int N, keytype* A)
{
	int k = 0;

		for (int i = 0; i < N; ++i) {

	    /* Invariant:
 			** - A[0:(k-1)] <= pivot; and
 			* - A[k:(i-1)] > pivot
		  */
		  const int ai = A[i];
	  	if (ai <= pivot) {
	    /* Swap A[i] and A[k] */
	      int ak = A[k];
	    	A[k++] = ai;
	    	 A[i] = ak;
	  	 }
		}
		return k;
}

int partition (keytype pivot, int N, keytype* A, keytype* B)
{
	  /* 0. if N < num of processors do serial partition*/
  int r = NUM_PROC;
  int c = N / NUM_PROC;
	int leftover = N % NUM_PROC;
  if(N < NUM_PROC)
  {
		paraPartition(pivot, N, A);
  }

	/* 2. do the partition */
	keytype lowPos[NUM_PROC];  //store k's :lowest idx of greater num
	#pragma omp parallel for private(i)   //[TODO: find if we need to share other vars]
	for(int i = 0; i < r; i++)
	{
		keytype count = i < leftover ? c+1:c ; //
		keytype start = i < leftover ? i * (c+1) : (leftover *(c+1) + (i - leftover) * c); // if leftover = 5 = i; count = c; start = 5*(c+1) 
		lowPos[i] = paraPartition(pivot, count, A+start);
	} 

	/* 3 count the boundary for each processor */  //[TODO: parallel it]
	int lowStart[NUM_PROC];
	int highStart[NUM_PROC];
	int lowsum = 0;
	#pragma omp parallel for
	for( int i = 0; i < NUM_PROC; i++)
	{
		lowStart[i] = lowsum;
		lowsum += lowPos[i];   // the final sum would be start of greater elements
	}
	int highsum = lowsum;
	#pragma parallel for
	for( int i = 0; i < NUM_PROC; i++)
	{
		highStart[i] = highsum;
		highsum += (i < leftover ? c+1 : c)  - lowPos[i];   // count number of greater elements of each processor's area.
	}

	/* 4. fill data into B*/
	#pragma omp parallel private(i)
	for(int i = 0; i < r ; i++)
	{
		
		for(int j = 0; j < lowPos[i] ; j++ )
		{
			B[j + lowStart[i]] = A[i*c + j];  
		}
		for(int j = 0; j < (i < leftover ? c+1 : c) - lowPos[i] ; j++)
		{
			keytype start = i < leftover ? i * (c+1) : (leftover *(c+1) + (i - leftover) * c);
			B[highStart[i] + j] = A[ start + lowPos[i] + j];  
		}
	}

	/* move data back to A*/
	#pragma omp parallel private(i)
	for ( int i = 0; i < r; i++)
	{
		int count = (i+1)*c > N ? (N - i*c) : c; 
		std::memcpy(A+i*c, B+i*c, count);
	}
  	return lowsum;  

}

void quickSort (int N, keytype* A, keytype* B)
{
  const int G = 1024; /* base :wqcase size, a tuning parameter */
  if (N < G)
    sequentialSort (N, A);
  else {
    // Choose pivot at random
    keytype pivot = A[rand () % N];

    // Partition around the pivot. Upon completion, n_less, n_equal,
    // and n_greater should each be the number of keys less than,
    // equal to, or greater than the pivot, respectively. Moreover, the array
    int n_le = partition (pivot, N, A, B);
		#pragma omp parallel sections
		{
			#pragma omp section
			{
    		quickSort (n_le, A, B);
			}
			#pragma omp section
			{
    		quickSort (N-n_le, A + n_le, B + n_le);
			}
		}
  }
}

void mySort (int N, keytype* A)
{
	keytype* B = (keytype*)malloc(N*sizeof(keytype));
  quickSort (N, A, B);
	free(B);
}

/* eof */
