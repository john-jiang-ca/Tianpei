/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include "cu_complex_operation.cuh"
#include <cuComplex.h>
#include <cuda.h>
//#include<helper_cuda.h>
//#include<helper_functions.h>
#include<time.h>
#include<cuda_profiler_api.h>
#include"common.h"




//pivot row normal1600ization
__global__ void normalizePivotRow( cuComplex *H, int index, int lda ) {
	int tid=threadIdx.x;
__shared__ cuComplex pivotValue;
if(tid==0)
{
//	printf("check the propagation matrix\n");
//	for(int count1=0;count1<N1;count1++)
//	{
//		for(int count2=0;count2<2*N1;count2++)
//		{
//			printf("%0.4f%+0.4fi ", H[IDC2D(count1,count2,lda)].x, H[IDC2D(count1,count2,lda)].y);
//		}
//		printf("\n");
//	}
	}
 if(tid<lda){

			if ( tid == 0 ) // First thread of each block loads pivotValue
			{
			pivotValue = H[ IDC2D( index, index, lda) ];
			}
			__syncthreads();

//			printf("the pivot value is %0.4f%+0.4fi :\n", pivotValue.x,pivotValue.y);
			H[ IDC2D( index, tid, lda )]=complex_div(H[ IDC2D( index, tid, lda )],pivotValue);
			__syncthreads();
 }

//printf("the thread Id of matrix iverse is:\n");
//printf("%d ", tid);
//printf("the row of the pivot H is:\n");
//printf("%0.4f%+0.4fi ", H[ IDC2D( index, tid, lda )].x,H[ IDC2D( index, tid, lda )].y );


}



//elements update
__global__ void linearMge( cuComplex *matrix, int index, int lda,int BLOCKNUM) {
	int tx=blockIdx.x*blockDim.x+threadIdx.x;
	int bid =threadIdx.x;
	int tid = threadIdx.y;
	 extern __shared__ cuComplex array[ ];
	__shared__ cuComplex zero;
		zero.x=0;
		zero.y=0;
//	extern __shared__ cuComplex matrixPivotValue[];
	cuComplex *matrixPivotValue=array;
	 cuComplex *multColumn=array+blockDim.x;
//			 int(int((lda)/2)/BLOCKNUM);
	 cuComplex *matrixRow=array+blockDim.x+lda;
if(tx<int(lda/2))
{
	if ( tx!=index ) {
		if(tid==0)
		{
	// Each block loads the value of the pivot Row to be substracted
	matrixPivotValue[bid] = matrix[ IDC2D( tx, index, lda )];
//	resultPivotValue = result[ IDC2D( index, x, lda )];
	matrix[ IDC2D(tx, index, lda )]=zero;
//	printf("the zeroing tx is %d:\n", tx);
		}
	}
	else
	{
		matrixPivotValue[bid]=zero;
	}
	__syncthreads();
	if(tid==0)
	{
//	printf("\n");
//	printf("the pivot column is:\n");
//	printf("%0.4f%+0.4fi, %d ",matrix[ IDC2D(tx, index, lda )].x,matrix[ IDC2D(tx, index, lda )].y, tx);
	}
	if(bid==0)
	{
	multColumn[ tid ] = matrix[ IDC2D( index, tid, lda )];
	}

	matrixRow[ IDC2D(bid,tid,lda) ] = matrix[ IDC2D( tx, tid, lda )];
//	resultRow[ ty ] = result[ IDC2D( y, x, lda )];
	__syncthreads();
//	newMatrixValue =matrix[ IDC2D( ty, x, lda )];
	if(tid!=index)
	{
matrix[ IDC2D(tx, tid, lda) ]=complex_sub(matrixRow[IDC2D(bid,tid,lda)],complex_mulcom( multColumn[tid],matrixPivotValue[bid]));
	}
	// Copy to the matrix
//	matrix[ IDC2D( ty, x, lda) ] = newMatrixValue;
	__syncthreads();

//	printf("the update value is:\n");
//		printf("%0.4f%+0.4fi ",matrix[ IDC2D( index, tid, lda )].x,matrix[ IDC2D( index, tid, lda )].y );
//	printf("the index of the whole matrix is %d:\n", tx);
//	printf("the index of the matrix in one block is: %d\n", bid);
}
}
__global__ void transfer(
		cuComplex *d_H,
		cuComplex *R_inv,
		int size
		)
{
	int bid=blockIdx.x*blockDim.x+threadIdx.x;
	int tid=threadIdx.y;
	if(bid>=0&&bid<size&&tid>=0&&tid<size)
	{
	R_inv[IDC2D(bid,tid,size)]=d_H[IDC2D(bid,(tid+size),2*size)];
	}
	__syncthreads();
//	if(bid==0&&tid==0)
//	{
//		printf("the result of transfer is:\n");
//		for(int count1=0;count1<size;count1++)
//		{//#define BLOCKNUM 4
	//Row switching
//			for(int count2=0;count2<size;count2++)
//			{
//				printf("%0.4f%+0.4fi ", R_inv[IDC2D(count1,count2,size)].x,R_inv[IDC2D(count1,count2,size)].y);
//			}
//			printf("\n");
//		}
//	}

}
__global__ void initial
(
		cuComplex *d_H,
		cuComplex *matrix,
		int size
		)
{
	int bid=blockIdx.x*blockDim.x+threadIdx.x;
	int tid=threadIdx.y;
	if(bid<size)
	{
	if(tid<size)
	{
	matrix[IDC2D(bid,tid,2*size)]=d_H[IDC2D(bid,tid,size)];
	}
	else if(tid==bid+size)
	{
		matrix[IDC2D(bid,tid,2*size)].x=1;
		matrix[IDC2D(bid,tid,2*size)].y=0;
	}
	else
	{
		matrix[IDC2D(bid,tid,2*size)].x=0;
		matrix[IDC2D(bid,tid,2*size)].y=0;
	}
	}
	__syncthreads();
//	if(bid==0&&tid==0)
//	{
//		printf("the result of initial is:\n");
//		for(int count1=0;count1<size;count1++)
//		{
//			for(int count2=0;count2<2*size;count2++)
//			{
//				printf("%0.4f%+0.4fi ", matrix[IDC2D(count1,count2,2*size)].x,matrix[IDC2D(count1,count2,2*size)].y);
//			}
//			printf("\n");
//		}
//	}
}
void MATRIX_INVERSE(
	cuComplex *H,  //input square matrix stored in row
	cuComplex *R,   //the inversion of the matrix H stored in row
	int row,        // the number of the rows
	int column      //the number of columns of the H_row
)
{
	int BLOCKNUM=16;
	cudaError_t error;
	int count1,count2;
	if(column<=8)
	{
		BLOCKNUM=1;
	}
	dim3 thread1(ceil(float(column)/float(BLOCKNUM)),2*column);
	dim3 thread2(ceil(float(column)/float(BLOCKNUM)),column);
	cuComplex *d_matrix;
	clock_t start, end;
	double duration;
	cudaMalloc((void**) &d_matrix, column*2*column*sizeof(cuComplex));
	start=clock();
	initial<<<BLOCKNUM,thread1>>>(H,d_matrix,column);
	end=clock();
	duration=double(end-start);
//	error=cudaDeviceSynchronize();
	if(error!=cudaSuccess)
	{
	printf("error=%s\n",cudaGetErrorString(cudaGetLastError()));
	}
	dim3 blockDim(ceil(float(column)/float(BLOCKNUM)),2*column);
	start=clock();
	for(count1=0; count1<column; count1++)
	{
   //shared memory to be changed
		normalizePivotRow<<<1,2*column>>>( d_matrix, count1, 2*column );

//		error=cudaDeviceSynchronize();
//        if(error!=cudaSuccess)
//        {
		printf("error=%s\n",cudaGetErrorString(cudaGetLastError()));
//        }
	linearMge<<<BLOCKNUM,blockDim,5000*sizeof(cuComplex)>>>( d_matrix, count1, 2*column,BLOCKNUM );
//	error=cudaDeviceSynchronize();
	if(error!=cudaSuccess)
//	{
	printf("error=%s\n",cudaGetErrorString(cudaGetLastError()));
//	}

	}
	end=clock();
	duration=double(end-start);
	start=clock();
  transfer<<<BLOCKNUM,thread2>>>(d_matrix,R,column);
	end=clock();
	duration=double(end-start);
  error=cudaDeviceSynchronize();
//	if(error!=cudaSuccess)
//	{
  printf("error=%s\n",cudaGetErrorString(cudaGetLastError()));
//	}
  cudaFree(d_matrix);
//	free(R);
//	free(H);

}




