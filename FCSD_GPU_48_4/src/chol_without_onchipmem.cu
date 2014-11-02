/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
//#include <assert.h>
//includes system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
//#include<cuda.h>
// includes CUDA
#include<cuda_runtime.h>
#include<cuda.h>
//#include<cutil.h>
//includes project
//#include<helper_cuda.h>
//#include<helper_functions.h>
#include"cu_complex_operation.cuh"
#include"common.h"
#include<cuComplex.h>
#include<time.h>
#include<cuda_profiler_api.h>
#define BLOCKNUM 32
 __global__ void chol_kernel_optimized
 (
		 cuComplex *R
 )
 {
//	 int MATRIX_SIZE=N1;
//	 extern __shared__ cuComplex R[];
      cuComplex pivot;
      cuComplex update;
	 int count1,count2,count3;
	 int tid=threadIdx.x;
//	 bid=blockIdx.x;
//	 tid=threadIdx.x;
	  cuComplex zero;
	 zero.x=0;
	 zero.y=0;
	 //cholesky factorization row by row
	 for (count1=0;count1<MATRIX_SIZE;count1++)
	 {
//		 R[tid]=matrix[IDC2D(count1,tid,MATRIX_SIZE)];
		 //pivoting step
			if(tid==0)
			{
				if(count1==0)
				{
					R[IDC2D(count1,count1,MATRIX_SIZE)].x=sqrt(R[IDC2D(count1,count1,MATRIX_SIZE)].x);
					R[IDC2D(count1,count1,MATRIX_SIZE)].y=0;
					printf("the first element is %0.4f%+0.4fi : \n",R[IDC2D(count1,count1,MATRIX_SIZE)].x,R[IDC2D(count1,count1,MATRIX_SIZE)].y );
				}
				else
				{
				pivot=zero;
				for(count2=0;count2<count1;count2++)
				{
			 pivot=complex_add(pivot,complex_mulcom(complex_conjugate(R[IDC2D(count2,count1,MATRIX_SIZE)]),R[IDC2D(count2,count1,MATRIX_SIZE)]));
				}
				 R[IDC2D(count1,count1,MATRIX_SIZE)]=complex_sub(R[IDC2D(count1,count1,MATRIX_SIZE)],pivot);
				 R[IDC2D(count1,count1,MATRIX_SIZE)].x=sqrt(R[IDC2D(count1,count1,MATRIX_SIZE)].x);
				 R[IDC2D(count1,count1,MATRIX_SIZE)].y=0;
//				 printf("the  %d diagonal elements is %0.4f%+0.4fi:\n",count1,R[IDC2D(count1,count1,MATRIX_SIZE)].x,R[IDC2D(count1,count1,MATRIX_SIZE)].y );
			}
			}
           __syncthreads();
// update the off-diagonal elements
				if(tid>count1)
				{
				  update=zero;
					for(count2=0;count2<count1;count2++)
					{
					update=complex_add(update,complex_mulcom(complex_conjugate(R[IDC2D(count2,count1,MATRIX_SIZE)]),R[IDC2D(count2,tid,MATRIX_SIZE)]));
					}
		R[IDC2D(count1,tid,MATRIX_SIZE)]=complex_sub(R[IDC2D(count1,tid,MATRIX_SIZE)],update);
		R[IDC2D(count1,tid,MATRIX_SIZE)]=complex_div(R[IDC2D(count1,tid,MATRIX_SIZE)],complex_conjugate(R[IDC2D(count1,count1,MATRIX_SIZE)]));
				}
				else if(tid<count1)
				{
					R[IDC2D(count1,tid,MATRIX_SIZE)]=zero;
				}
		      __syncthreads();
				printf("the %d %d element is %0.4f%+0.4fi ", count1,tid, R[IDC2D(count1,tid,MATRIX_SIZE)].x,R[IDC2D(count1,tid,MATRIX_SIZE)].y);
				printf("\n");
			}

//	 if(tid>=0&&tid<MATRIX_SIZE)
//	 {
//		 for(count1=0;count1<MATRIX_SIZE;count1++)
//		 {
//		 matrix[IDC2D(count1,tid,MATRIX_SIZE)]=R[IDC2D(count1,tid,MATRIX_SIZE)];
//		 }
//	 }
//	 __syncthreads();

 }


void chol_without_onchip(cuComplex *d_U
)
{


	//int MATRIX_SIZE=sizeof(U)/;
	//int m=U->size2;
//	int count;
	int count1;
	int count2;
//int MATRIX_SIZE=N1;
//allocate computation space
	cudaError_t error;
//	cuComplex *d_U;
//	error=cudaMalloc((void**) &d_U, MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
//	if(error !=cudaSuccess){
//		printf("cudaMalloc d_U returned error code %d, line(%d)\n", error, __LINE__);
//	}
	//data transmission from CPU to GPU


	//data transmission from CPU to GPU


//	int sharedMem=MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex);
	//dim3 BlockId;
	clock_t start, end;
	start=clock();
		int threadNum=MATRIX_SIZE;
		int blockNum=1;
		int sharedMem=2000*sizeof(cuComplex);
//		int *threadID;
//		int *d_threadID;
		cudaProfilerStart();
//		threadID=(int*)malloc(sizeof(int)*threadNum);
//		cudaMalloc((void**)&d_threadID, threadNum*sizeof(int));
//		cudaMemcpy(d_threadID,threadID,sizeof(int)*threadNum,cudaMemcpyHostToDevice);
//		error=cudaMemcpy(d_U, U, sizeof(cuComplex)*MATRIX_SIZE*MATRIX_SIZE, cudaMemcpyHostToDevice );
//		if(error!=cudaSuccess){
//			printf("cudaMemcpy U to d_U returned error code %d, line(%d)\n", error, __LINE__);
//		}

		chol_kernel_optimized<<<blockNum,threadNum,sharedMem>>>(d_U);
//		cudaError_t error;
		error=cudaDeviceSynchronize();
		if(error!=cudaSuccess)
		{
		printf("%s\n", cudaGetErrorString(cudaGetLastError()));
		}
//		cudaMemcpy(threadID, d_threadID, sizeof(int)*threadNum, cudaMemcpyDeviceToHost );
//		for(count1=0;count1<threadNum; count1++)
//		{
//		printf("%d", threadID[count1]);
//		printf("\n");
//		}
//		printf("%d",sizeof(cuComplex));
//		error=cudaMemcpy(pR,d_U,sizeof(cuComplex)*MATRIX_SIZE*MATRIX_SIZE, cudaMemcpyDeviceToHost );
//		if(error!=cudaSuccess)
//		{
//			printf("cudaMemcpy d_U to U returned error code %d, line(%d)\n", error, __LINE__);
//		}
		cudaProfilerStop();
//	}
	end=clock();
//	 *durationD=(double)(end-start)/CLOCKS_PER_SEC;
	//data transmission from GPU to CPU
//	for(count1=0; count1<MATRIX_SIZE; count1++)
//	{
//		for(count2=0; count2<MATRIX_SIZE; count2++)
//		{
//			pR[count1*MATRIX_SIZE+count2]=U[count1*MATRIX_SIZE+count2];
//
//
//		}
//	}
//	U[0].x=*m;
//	pR[0].y=*m;
//    free(y);
//    free(m);
//    free(threadID);
//    cudaFree(d_threadID);

}

