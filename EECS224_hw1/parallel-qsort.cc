/**
 *  \file parallel-qsort.cc
 *
 *  \brief Implement your parallel quicksort algorithm in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

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
int partition (keytype pivot, int N, keytype* A)
{
  int k = 0;

	  /* 0. if N < num of processors do serial partition*/
  int r = NUM_PROC;
  int c = N / NUM_PROC + ((N%NUM_PROC==0) ? 0 : 1);
  if(N < NUM_PROC)
  {
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

	/* 1. create 2 arrays */
	keytype **lower = (keytype **)malloc(r * sizeof(keytype*));
	keytype **higher = (keytype **)malloc(r * sizeof(keytype*));
	for ( int i=0; i < r; i ++)
	{
		lower[i] = (keytype*)malloc((c+1)*sizeof(keytype));  //add 1 more square infront of array to store count
		higher[i] = (keytype*)malloc((c+1)*sizeof(keytype));
	}

	/* 2. do the partition */
	#pragma omp parallel for shared(lower, higher) private(i)   //[TODO: find if we need to share other vars]
	for(int i = 0; i < r; i++)
	{
		int lowIndex = 1;
		int highIndex = 1;
		for (int j = i * c; j < (i+1)*c && j < N; j++)
		{
			if(A[j] < pivot)
				lower[i][lowIndex++] = A[j];
			else if(A[j] > pivot)
				higher[i][highIndex++] = A[j];
		}

		lower[i][0] = lowIndex - 1;
		higher[i][0] = highIndex - 1;
	}

	/* 3 count the boundary for each processor */  //[TODO: parallel it]
	keytype* lowStart = (keytype*) malloc(r * sizeof(keytype));
	keytype* highTail = (keytype*) malloc(r*sizeof(keytype));
	keytype curLow = 0;
	keytype curHigh = N-1;
	for(int i = 0 ; i < r; i++)
	{
		lowStart[i] = curLow;
		highTail[r - i - 1] = curHigh;

		curLow += lower[i][0];
		curHigh -= higher[r - i - 1][0];
	}

	/* 4. fill back A*/
	int subgap = (curHigh - curLow + 1) / r + ((curHigh - curLow + 1)%r == 0? 0:1);

	#pragma omp parallel shared(lowStart, highTail) private(i)
	for(int i = 0; i < r ; i++)
	{
		
		for(int j = 0; j < lower[i][0] ; j++ )
		{
			A[j + lowStart[i]] = lower[i][j + 1];  //skip lower[i][0]
		}
		for(int j = 0; j < higher[i][0]; j++)
		{
			A[highTail[i] - j] = higher[i][higher[i][0] - j];  // omit higher[i][0]
		}

		/* gap */
		if(curHigh >= curLow)
		{
			for(int j = i*subgap + curLow ; j < (i+1)*subgap + curLow && j <= curHigh; j++)
				A[j] = pivot;
		}
	}

	/* free the spaces */
	for ( int i=0; i < r; i ++)
	{
		free(lower[i]); 
		free(higher[i]);
	}
	
	free(lower);
	free(higher);
	free(lowStart);
	free(highTail);


  	return curHigh + 1;  

}

void quickSort (int N, keytype* A)
{
  const int G = 1024; /* base case size, a tuning parameter */
  if (N < G)
    sequentialSort (N, A);
  else {
    // Choose pivot at random
    keytype pivot = A[rand () % N];

    // Partition around the pivot. Upon completion, n_less, n_equal,
    // and n_greater should each be the number of keys less than,
    // equal to, or greater than the pivot, respectively. Moreover, the array
    int n_le = partition (pivot, N, A);
		#pragma omp parallel sections
		{
			#pragma omp section
			{
    		quickSort (n_le, A);
			}
			#pragma omp section
			{
    		quickSort (N-n_le, A + n_le);
			}
		}
  }
}

void mySort (int N, keytype* A)
{
  quickSort (N, A);
}

/* eof */
