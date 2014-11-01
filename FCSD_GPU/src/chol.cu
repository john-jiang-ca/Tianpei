/*
 /*
 *  THIS FUNCTION INPLEMENT CHOLESKY FACTORIZATION
 *
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
//#define BLOCKNUM 32
 __global__ void chol_kernel_optimized
 (
		 cuComplex *R,
//		 cuComplex *pivot,   //using for the data communication among all the blocks
		 int index
 )
 {
//	 int MATRIX_SIZE=N1;
	 int count1,count2;
	 int tid=threadIdx.x;
      __shared__ cuComplex update;
      extern __shared__ cuComplex array[];
	 cuComplex *matrixColumn=array;
	 __shared__  cuComplex zero;
	 zero.x=0;
	 zero.y=0;
	 //cholesky factorization row by row
//	 for (count1=0;count1<MATRIX_SIZE;count1++)
//	 {
//		 R[tid]=matrix[IDC2D(count1,tid,MATRIX_SIZE)];
//transfer the data into share array in one block
	 if(tid>=0&&tid<(MATRIX_SIZE-index))
	 {
	 for(count1=0;count1<=index;count1++)
	 {
          matrixColumn[IDC2D(count1,tid,(MATRIX_SIZE-index))]=R[IDC2D(count1,(tid+index),MATRIX_SIZE)];


	 }
//	 printf("this is working!!");
	 }
 __syncthreads();
// printf("%d ", tid);
// if(tid==0)
// {
//	 printf("hey I am here!!");
//	 for(count1=0;count1<=index;count1++)
//	 {
//		 for(count2=0;count2<MATRIX_SIZE;count2++)
//		 {
//          printf("%0.4f%+0.4fi:",R[IDC2D(count1,count2,(MATRIX_SIZE))].x,R[IDC2D(count1,count2,(MATRIX_SIZE))].y );
//		 }
//		 printf("\n");
//	 }
//}
// __syncthreads();
 //pivoting step
			if(tid==0)
			{
				if(index==0)
				{

					matrixColumn[IDC2D(0,0,(MATRIX_SIZE-index))].x=sqrt(matrixColumn[IDC2D(0,0,(MATRIX_SIZE-index))].x);
					matrixColumn[IDC2D(0,0,(MATRIX_SIZE-index))].y=0;
					R[IDC2D(0,0,MATRIX_SIZE)]=matrixColumn[IDC2D(0,0,(MATRIX_SIZE-index))];
//					printf("the first element is %0.4f%+0.4fi : \n",matrixColumn[IDC2D(index,index,MATRIX_SIZE)].x,matrixColumn[IDC2D(index,index,MATRIX_SIZE)].y );

				}
				else
				{
				update=zero;
//the addition of the column per block
				for(count2=0;count2<index;count2++)
				{
			 update=complex_add(update,complex_mulcom(complex_conjugate(matrixColumn[IDC2D(count2,0,(MATRIX_SIZE-index))]),matrixColumn[IDC2D(count2,0,(MATRIX_SIZE-index))]));
				}
				matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))]=complex_sub(matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))],update);
				matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))].x=sqrt(matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))].x);
				matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))].y=0;
				R[IDC2D(index,index,MATRIX_SIZE)]=matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))];
//				pivot[blockIdx.x]=update;
//				 printf("the  %d diagonal elements is %0.4f%+0.4fi:\n",count1,R[IDC2D(count1,count1,MATRIX_SIZE)].x,R[IDC2D(count1,count1,MATRIX_SIZE)].y );
			}
			}
__syncthreads();
// update the off-diagonal elements
if(tid>0&&tid<(MATRIX_SIZE-index))
{
if(index==0)
{
	R[IDC2D(index,(tid+index),MATRIX_SIZE)]=complex_div(matrixColumn[IDC2D(index,tid,(MATRIX_SIZE-index))],complex_conjugate(matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))]));
}
else
{
for(count2=0;count2<index;count2++)
{
matrixColumn[IDC2D(index,tid,(MATRIX_SIZE-index))]=complex_sub(matrixColumn[IDC2D(index,tid,(MATRIX_SIZE-index))],complex_mulcom(complex_conjugate(matrixColumn[IDC2D(count2,0,(MATRIX_SIZE-index))]),matrixColumn[IDC2D(count2,tid,(MATRIX_SIZE-index))]));

}
R[IDC2D(index,(tid+index),MATRIX_SIZE)]=complex_div(matrixColumn[IDC2D(index,tid,(MATRIX_SIZE-index))],complex_conjugate(matrixColumn[IDC2D(index,0,(MATRIX_SIZE-index))]));

//printf("hey attention!:\n");
//printf("%0.4f%+0.4fi ",matrixColumn[IDC2D(index,tid,(pitch))].x,matrixColumn[IDC2D(index,tid,(MATRIX_SIZE-index))].y );
}
}
__syncthreads();


//				printf("the %d %d element is %0.4f%+0.4fi ", count1,tid, R[IDC2D(count1,tid,MATRIX_SIZE)].x,R[IDC2D(count1,tid,MATRIX_SIZE)].y);
//printf("I am working !!\n");
 }

__global__ void zeroing(
		cuComplex *H
		)
{
//	int MATRIX_SIZE=N1;
 int tid=threadIdx.x;
 int count1,count2;
 for(count1=0;count1<MATRIX_SIZE;count1++)
 {
 if(tid>=0&&tid<=(count1-1))
 {
	 H[IDC2D(count1,tid,(MATRIX_SIZE))].x=0;
	 H[IDC2D(count1,tid,(MATRIX_SIZE))].y=0;
 }
 }
// if(tid==0)
// {
// for(count1=0;count1<MATRIX_SIZE;count1++)
// {
//	 for(count2=0;count2<MATRIX_SIZE;count2++)
//	 {
//		 printf("%0.4f%+0.4fi ", H[IDC2D(count1,count2,(MATRIX_SIZE))].x, H[IDC2D(count1,count2,(MATRIX_SIZE))].y);
//	 }
//	 printf("\n");
// }
// }
}
void chol(cuComplex *d_U

)
{


//	int MATRIX_SIZE=N1;
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
	//dim3 BlockId;printf("hey attention!:\n");
//	clock_t start, end;
//	start=clock();

//		int *threadID;
//		int *d_threadID;
//		cudaProfilerStart();
//		threadID=(int*)malloc(sizeof(int)*threadNum);
//		cudaMalloc((void**)&d_threadID, threadNum*sizeof(int));
//		cudaMemcpy(d_threadID,threadID,sizeof(int)*threadNum,cudaMemcpyHostToDevice);
//		error=cudaMemcpy(d_U, U, sizeof(cuComplex)*MATRIX_SIZE*MATRIX_SIZE, cudaMemcpyHostToDevice );
//		if(error!=cudaSuccess){
//			printf("cudaMemcpy U to d_U returned error code %d, line(%d)\n", error, __LINE__);
//		}
//	   double duration;
//	   clock_t start, end;
//	   start=clock();
   for(count1=0;count1<MATRIX_SIZE;count1++)
   {
		int threadNum=(MATRIX_SIZE-count1);
		int blockNum=1;
		int sharedMem=4000*sizeof(cuComplex);
		chol_kernel_optimized<<<blockNum,threadNum,sharedMem>>>(d_U,count1);
//		error=cudaDeviceSynchronize();
				if(error!=cudaSuccess)
				{
				printf("%s\n", cudaGetErrorString(cudaGetLastError()));
				}
   }
//end=clock();
//duration=double(end-start);
//printf("%0.4f ", duration);
//printf("\n");
//  start=clock();
   zeroing<<<1,MATRIX_SIZE>>>(d_U);
	error=cudaDeviceSynchronize();
//	   end=clock();
//	   duration=double(end-start);
//	   printf("%0.4f ", duration);
//	   printf("\n");
			if(error!=cudaSuccess)
			{
			printf("%s\n", cudaGetErrorString(cudaGetLastError()));
			}

//		cudaError_t error;

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
//		cudaProfilerStop();
//	}
//	end=clock();
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
