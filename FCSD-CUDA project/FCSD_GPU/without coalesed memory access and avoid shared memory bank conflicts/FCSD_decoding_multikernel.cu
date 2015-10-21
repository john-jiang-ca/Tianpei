/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
/*
 * this function implement fixed complexity sphere decoding
 * INPUT:
 * y: received signal
 * H: permuted propagation matrix
 * M: modulation scheme, (2: BPSK 4: QPSK, 16: 16QAM, 64: 64QAM)
 * psymbolconstellation: the symbol constellation
 * OUTPUT:
 * s: detection result
 * Eu: Euclidean distance
 */
#include <stdio.h>
#include <stdlib.h>
#include<cuComplex.h>
#include <string.h>
#include<math.h>
#include<cuda_runtime.h>
#include<cuda.h>
#include<cublas_v2.h>
#include<time.h>
#include"cu_complex_operation.cuh"
#include"common.h"
#include<cuda_runtime.h>
#include<cuda.h>
#include<cudaProfiler.h>
#include<cuda_profiler_api.h>
#include<device_functions.h>
#define threadNum 1024
#define blockNum 1
#define stride 1
/*
 * in this version I applied colesced memory accesss to all the vector and matrix, with all the matrix stored in row major different threads reading one column
 * with all the matrix stored in column major, different threads reading one row
 * uitilize consecutive computation power
 */
__constant__ cuComplex d_R[MATRIX_SIZE*MATRIX_SIZE],d_psymbolconstellation[16],d_constant_shat[MATRIX_SIZE];
__constant__ int d_list[MATRIX_SIZE];
__global__ void FEpath(
//		cuComplex *d_R,  //upper triangular matrix after cholesky factorization
//		cuComplex *d_constant_shat,  //unconstrained estimation of transmitted symbol vector s
//		cuComplex *s_matrix_share,
		cuComplex *s_potential_matrix,   //the matrix use to store all the solution candidates from all the blocks
		int *s_sub_index,   //full factorial index matrix
		int rho,
		int pathNum,
//		cuComplex *s,  //decoding results
		float *Eu,  //Euclidean distance
//		int pitch_R,    //the number of transmit antennas
//		int pitch_index,    //the number of receive antennas
//		int pitch_p,
		int M,    //modulation scheme
//		int threadNum,    //number of threads
//		int *d_list,     //the permutation list
//		cuComplex *d_psymbolconstellation,
		int index





		)
{
	//need to consider the resource allocation
	int tx=blockIdx.x*blockDim.x+threadIdx.x;     //if the path number is small we can allocate the kernel into one block so that we can use the shared memory
int tid=threadIdx.x;
int Nt=MATRIX_SIZE;

//allocate shared memory
	extern __shared__ cuComplex array[];

	error_t error;
	int count1, count2,count3,count4;
	__shared__ float d;    //the minimum distance unit between the signal constellation, the distance is usually 2d
	__shared__ cuComplex alpha, beta;
	alpha.x=1;alpha.y=0; beta.x=0; beta.y=0;
	__shared__ cuComplex s_hat_share[MATRIX_SIZE];
	__shared__ cuComplex R_share[MATRIX_SIZE*MATRIX_SIZE];
//    __shared__ cuComplex s_temp[threadNum];
	cuComplex s_temp;
	cuComplex Eu_norm_share;
//    __shared__ cuComplex Eu_norm_share[threadNum];
    cuComplex *R_Eu_share=(cuComplex*)malloc(Nt*sizeof(cuComplex));
    if(tid>=0&&tid<MATRIX_SIZE)
    {

		for(count1=0;count1<MATRIX_SIZE;count1++)
		{
			R_share[IDC2D(tid,count1,MATRIX_SIZE)]=d_R[IDC2D(tid,count1,MATRIX_SIZE)];
		}


		s_hat_share[tid]=d_constant_shat[tid];
    }

	__syncthreads();

	Eu[blockNum*threadNum*index+tx]=0;
for (count1=Nt-1; count1>=0; count1--)
{

		if (count1<Nt-rho)
		{
			s_temp=s_hat_share[count1];
            #pragma unroll
			for (count2=count1+1;count2<Nt; count2++)
			{
				s_temp=complex_add(s_temp,complex_mulcom(complex_div(R_share[IDC2D(count1,count2,MATRIX_SIZE)],R_share[IDC2D(count1,count1,MATRIX_SIZE)]),(complex_sub(s_hat_share[count2],s_potential_matrix[IDC2D(count2,(index*blockNum*threadNum+tx),pathNum/4)]))));
			}
	if(M==2)   //BPSK
	{

    	 d=sqrt(float(float(1)/float(Nt)));
    	 if(s_temp.x>0)
    	 {
    		  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=d;
    		  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=0;
    	 }
    	 else
    	 {
    		  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=(-d);
    		  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=0;
    	 }
	}
	else if(M==4)   //QPSK
	{
		int *result=(int*)malloc(sizeof(int));

		memset(result,0,sizeof(int));
		float *distance=(float*)malloc(M*sizeof(float));
		memset(distance,0,M*sizeof(int));
	    	d=sqrt(float(float(1)/float(Nt)));
           for(count2=0;count2<M;count2++)
           {
        	   switch(count2)
        	   {
        	   case 0:
        	   distance[count2]=sqrt(powf(s_temp.x-(-d),2)+pow(s_temp.y-0,2)); break;
        	   case 1:
        	   distance[count2]=sqrt(powf(s_temp.x-0,2)+pow(s_temp.y-(-d),2)); break;
        	   case 2:
        	   distance[count2]=sqrt(pow(s_temp.x-(d),2)+pow(s_temp.y-0,2)); break;
        	   case 3:
        	   distance[count2]=sqrt(pow(s_temp.x-0,2)+pow(s_temp.y-d,2)); break;
               default:
                   #if __CUDA_ARCH__ >=300
            	   printf("result error code %d\n", error);
                   #endif
            	   break;
        	   }
           }
           float  mini_distance;
  	     int mini_index;
  	     mini_distance=distance[0];
  	     mini_index=1;
  	     for(count3=0;count3<M;count3++)
  	     {
  	    	if(distance[count3]<mini_distance)
  	    	{
  	    		mini_distance=distance[count3];
  	    		mini_index=count3+1;
  	    	}
  	     }

           switch (mini_index)
           {
           case 1:
        	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=-d;
        	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=0;
        	break;
           case 2:
           	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=0;
           	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=-d;
           	break;
           case 3:
           	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=d;
           	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=0;
           	break;
           case 4:
       	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=0;
       	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=d;
       	  break;
           default:
               #if __CUDA_ARCH__ >=300
        	   printf("result error code %d\n", error);
               #endif
           }
           free(distance);
           free(result);
	}
	else if(M==16)  //16QAM
	{

	    	 d=sqrt(float(3)/(2* (float)(Nt*(M-1))));
	     if(s_temp.x<(-2*d))
	    	{
	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=(-3*d);
	    	}
	    else if(s_temp.x>(2*d))
	    	{
	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=(3*d);
	    	}
	    else if(s_temp.x>=0&&s_temp.x<=2*d)
	    	{
	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=d;
	    	}
	    else if(s_temp.x>=(-2*d)&&s_temp.x<=0)
	        {
	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].x=(-d);
	    	}




	    if(s_temp.y<(-2*d))
	    	{
	    	    s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=(-3*d);
	    	}
	   else if(s_temp.y>(2*d))
	    	 {
	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=(3*d);
	         }
	   else if(s_temp.y>=0&&s_temp.y<=(2*d))
	        {
	    	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=d;
	    	}
	    else if(s_temp.y>=(-2*d)&&s_temp.y<=0)
	       {
	         s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)].y=(-d);
	       }
	}
	else if(M==64)   //64QAM
		{

	    	 d=sqrt(3/(2* (float)(Nt*(M-1))));
	     }
		}
	else
	{
		 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)]=d_psymbolconstellation[s_sub_index[IDC2D((Nt-count1-1),(index*blockNum*threadNum+tx),pathNum/4)]];

	}
		R_Eu_share[count1]=complex_sub( s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pathNum/4)],s_hat_share[count1]);
		Eu_norm_share=beta;
        #pragma unroll
		for(count3=count1;count3<MATRIX_SIZE;count3++)
		{
		Eu_norm_share=complex_add(Eu_norm_share,complex_mulcom(R_share[IDC2D(count1,count3,MATRIX_SIZE)],R_Eu_share[count3]));
		}
		Eu[(index*blockNum*threadNum+tx)]=Eu[(index*blockNum*threadNum+tx)]+pow(Eu_norm_share.x,2)+pow(Eu_norm_share.y,2);

	}
__syncthreads();
free(R_Eu_share);


}

//host

void FCSD_decoding(
		cuComplex *R,  //upper triangular matrix after cholesky factorization store in device side
//		cuComplex *s_sub, //the sub brute force rho vector matrix
		cuComplex *d_s_hat,  //unconstrained estimation of transmitted symbol vector s
		cuComplex *s_kernel,  //quantization of estimation ,decoding results
//		cuComplex *Eu,  //Euclidean distance
		int Nt,    //the number of transmit antennas
		int Nr,    //the number of receive antennas
		int M,    //modulation scheme
		int *list,   //the permutation list
		cuComplex *psymbolconstellation //the symbol constellation
		)
{
//brute force search determine the vector results of the full expansion
	int rho=ceil(sqrt(Nt)-1);
	int count1,count2;
//	cuComplex *ss;
//	ss=(cuComplex*)malloc(MATRIX_SIZE*sizeof(cuComplex));
//	cuComplex *s_sub;
//	s_sub=(cuComplex*)malloc(pow(M,rho)*rho*sizeof(cuComplex));   //all the possible full expansion sub vector
	int  pathNum;
	pathNum=pow(M,rho);
//	int *d_s_sub_index;
	int *s_sub_index=(int*)calloc(1,rho*pow(M,rho)*sizeof(int));
	fullfact(rho,M,s_sub_index);    //get  the indexes of all the possible rho length symbol vectors
//	int blockNum=BLOCK_NUM;   //determined by the path number
//	int pathNum=pow(M,rho);  //number of search path
//	int threadNum=ceil(pathNum/(blockNum*stride)); //determined by the path number
    float *Eu;
	Eu=(float*)calloc(1,blockNum*sizeof(float));
//	cuComplex *s_potential_matrix=(cuComplex*)calloc(1,pathNum*Nt*sizeof(cuComplex));
	cuComplex *s_hat=(cuComplex*)calloc(1,Nt*sizeof(cuComplex));
//	cuComplex  *d_s_potential_matrix;
	cudaMemcpy(s_hat,d_s_hat,Nt*sizeof(MATRIX_SIZE),cudaMemcpyDeviceToHost);

//	int *j;
//	j=(int*)malloc(sizeof(int));
	cublasHandle_t handle;
		cublasStatus_t ret;
		cudaError_t error;
		size_t pitch_R,pitch_potential,pitch_index;
		ret=cublasCreate(&handle);
//	    error=cudaMalloc((void**) &d_R, MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
//	    error=cudaMalloc((void**) &d_s_sub_index, rho*pathNum*sizeof(int));
//		error=cudaMalloc((void**) &d_s_potential_matrix, pathNum*Nt*sizeof(cuComplex));
//		error=cudaMalloc((void**) &d_list, MATRIX_SIZE*sizeof(int));
//		error=cudaMalloc((void**) &d_Eu, pathNum*sizeof(float));
//        error=cudaMalloc((void**) &d_psymbolconstellation, M*sizeof(cuComplex));
        clock_t start, end;
        start=clock();
        cuComplex *R_constant=(cuComplex*)calloc(1,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
    	cudaMemcpyToSymbol(d_R, R, MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),0,cudaMemcpyHostToDevice);
    	cudaMemcpyFromSymbol(R_constant,d_R,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),0,cudaMemcpyDeviceToHost);
    	printf("cao ni ma!!\n");
    	for(count1=0;count1<MATRIX_SIZE;count1++)
    	{
    		for(int count2=0;count2<MATRIX_SIZE;count2++)
    		{
    			printf("%0.4f%+0.4fi ", R_constant[IDC2D(count1,count2,MATRIX_SIZE)].x,R_constant[IDC2D(count1,count2,MATRIX_SIZE)].y);
    		}
    		printf("\n");
    	}
    	cudaMemcpyToSymbol(d_psymbolconstellation, psymbolconstellation, M*sizeof(cuComplex),0,cudaMemcpyHostToDevice);
    	cudaMemcpyToSymbol(d_constant_shat, s_hat, MATRIX_SIZE*sizeof(cuComplex),0,cudaMemcpyHostToDevice);
    	cudaMemcpyToSymbol(d_list, list, MATRIX_SIZE*sizeof(int),0,cudaMemcpyHostToDevice);
//        error=cudaMemcpy(d_R,R,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyDeviceToDevice);
//		error=cudaMemcpy(d_psymbolconstellation, psymbolconstellation, M*sizeof(cuComplex),cudaMemcpyHostToDevice);
//		error=cudaMemcpy(d_s_sub_index, s_sub_index,rho*pathNum*sizeof(int),cudaMemcpyHostToDevice);
//		error=cudaMemcpy(d_list, list, Nt*(sizeof(int)),cudaMemcpyHostToDevice);
	 int sharedMem;
    sharedMem=1*sizeof(cuComplex);
   float duration;
   cuComplex *s_potential1,*s_potential2,*s_potential3,*s_potential4;
   cuComplex *d_s_potential1,*d_s_potential2,*d_s_potential3,*d_s_potential4;
	float *Eu1,*Eu2,*Eu3,*Eu4;
	float *d_Eu1,*d_Eu2,*d_Eu3,*d_Eu4;
	int *sub_index1,*sub_index2,*sub_index3,*sub_index4;
	int *d_sub_index1,*d_sub_index2,*d_sub_index3,*d_sub_index4;
    cudaHostAlloc((void**) &s_potential1,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex), cudaHostAllocDefault);
    cudaHostAlloc((void**) &s_potential2,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex), cudaHostAllocDefault);
    cudaHostAlloc((void**) &s_potential3,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex), cudaHostAllocDefault);
    cudaHostAlloc((void**) &s_potential4,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Eu1,(pathNum/4)*sizeof(float), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Eu2,(pathNum/4)*sizeof(float), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Eu3,(pathNum/4)*sizeof(float), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Eu4,(pathNum/4)*sizeof(float), cudaHostAllocDefault);
    cudaHostAlloc((void**) &sub_index1,(pathNum/4)*sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &sub_index2,(pathNum/4)*sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &sub_index3,(pathNum/4)*sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &sub_index4,(pathNum/4)*sizeof(int), cudaHostAllocDefault);
    cudaMalloc((void**) &d_Eu1,(pathNum/4)*sizeof(float));
    cudaMalloc((void**) &d_Eu2,(pathNum/4)*sizeof(float));
    cudaMalloc((void**) &d_Eu3,(pathNum/4)*sizeof(float));
    cudaMalloc((void**) &d_Eu4,(pathNum/4)*sizeof(float));
    cudaMalloc((void**) &d_s_potential1,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex));
    cudaMalloc((void**) &d_s_potential2,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex));
    cudaMalloc((void**) &d_s_potential3,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex));
    cudaMalloc((void**) &d_s_potential4,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex));
    cudaMalloc((void**) &d_sub_index1,rho*(pathNum/4)*sizeof(int));
    cudaMalloc((void**) &d_sub_index2,rho*(pathNum/4)*sizeof(int));
    cudaMalloc((void**) &d_sub_index3,rho*(pathNum/4)*sizeof(int));
    cudaMalloc((void**) &d_sub_index4,rho*(pathNum/4)*sizeof(int));
    for(count1=0;count1<int(pathNum/4);count1++)
    {
    	for( count2=0;count2<rho;count2++)
    	{
         sub_index1[IDC2D(count2,count1,pathNum/4)]= s_sub_index[IDC2D(count2,count1,pathNum)];
    	}
    }
    for(count1=0;count1<int(pathNum/4);count1++)
    {
    	for(count2=0;count2<rho;count2++)
    	{
         sub_index2[IDC2D(count2,count1,pathNum/4)]= s_sub_index[IDC2D(count2,count1+int(pathNum/4),pathNum)];
    	}
    }
    for(count1=0;count1<int(pathNum/4);count1++)
    {
    	for(count2=0;count2<rho;count2++)
    	{
         sub_index3[IDC2D(count2,count1,pathNum/4)]= s_sub_index[IDC2D(count2,count1+2*int(pathNum/4),pathNum)];
    	}
    }
    for(count1=0;count1<int(pathNum/4);count1++)
    {
    	for(count2=0;count2<rho;count2++)
    	{
         sub_index4[IDC2D(count2,count1,pathNum/4)]= s_sub_index[IDC2D(count2,count1+3*int(pathNum/4),pathNum)];
    	}
    }

    cudaStream_t stream1,stream2,stream3,stream4;
    cudaStreamCreate(&stream1);cudaStreamCreate(&stream2);cudaStreamCreate(&stream3);cudaStreamCreate(&stream4);
   for(count1=0;count1<stride;count1++)
   {
	cudaMemcpyAsync(d_sub_index1,sub_index1,rho*(pathNum/4)*sizeof(int),cudaMemcpyHostToDevice,stream1);
	FEpath<<<blockNum, threadNum,sharedMem,stream1>>>(d_s_potential1,d_sub_index1,rho,pathNum, d_Eu1, M,count1);
	error=cudaDeviceSynchronize();
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpyAsync(s_potential1,d_s_potential1,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex),cudaMemcpyDeviceToHost,stream1);
	cudaMemcpyAsync(Eu1,d_Eu1,(pathNum/4)*sizeof(float),cudaMemcpyDeviceToHost,stream1);

	cudaMemcpyAsync(d_sub_index2,sub_index2,rho*(pathNum/4)*sizeof(int),cudaMemcpyHostToDevice,stream2);
	FEpath<<<blockNum, threadNum,sharedMem,stream2>>>(d_s_potential2,d_sub_index2,rho,pathNum, d_Eu2, M,count1);
	error=cudaDeviceSynchronize();
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpyAsync(s_potential2,d_s_potential2,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex),cudaMemcpyDeviceToHost,stream2);
	cudaMemcpyAsync(Eu2,d_Eu2,(pathNum/4)*sizeof(float),cudaMemcpyDeviceToHost,stream2);

	cudaMemcpyAsync(d_sub_index3,sub_index3,rho*(pathNum/4)*sizeof(int),cudaMemcpyHostToDevice,stream3);
	FEpath<<<blockNum, threadNum,sharedMem,stream3>>>(d_s_potential3,d_sub_index3,rho,pathNum, d_Eu3, M,count1);
	error=cudaDeviceSynchronize();
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpyAsync(s_potential3,d_s_potential3,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex),cudaMemcpyDeviceToHost,stream3);
	cudaMemcpyAsync(Eu3,d_Eu3,(pathNum/4)*sizeof(float),cudaMemcpyDeviceToHost,stream3);

	cudaMemcpyAsync(d_sub_index4,sub_index4,rho*(pathNum/4)*sizeof(int),cudaMemcpyHostToDevice,stream4);
	FEpath<<<blockNum, threadNum,sharedMem,stream4>>>(d_s_potential4,d_sub_index4,rho,pathNum, d_Eu4, M,count1);
	cudaMemcpyAsync(s_potential4,d_s_potential4,MATRIX_SIZE*(pathNum/4)*sizeof(cuComplex),cudaMemcpyDeviceToHost,stream4);
	cudaMemcpyAsync(Eu4,d_Eu4,(pathNum/4)*sizeof(float),cudaMemcpyDeviceToHost,stream4);
	error=cudaDeviceSynchronize();
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
   }
//	error=cudaMemcpy(s_potential_matrix,d_s_potential_matrix, Nt*sizeof(cuComplex)*pathNum,cudaMemcpyDeviceToHost);
	end=clock();
//	cudaProfilerStop();
	duration=double(end-start);

	printf("hey %0.4f ", duration);
	printf("\n");
//    memcpy(Eu+0,Eu1,pathNum/4*sizeof(float));
//    memcpy(Eu+pathNum/4,Eu2,pathNum/4*sizeof(float));
//    memcpy(Eu+pathNum/2,Eu3,pathNum/4*sizeof(float));
//    memcpy(Eu+(pathNum*3)/4,Eu4,pathNum/4*sizeof(float));

//	printf("Eu_num is %d", Eu_num);
//    error=cudaMemcpy(s_potential_matrix,d_s_potential_matrix,pathNum*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);

//    printf("all the potential symbol vector is:\n");
//    for(count1=0;count1<pathNum;count1++)
//    {
//    	for(int count2=0;count2<Nt;count2++)
//    	{
//    		printf("%0.4f%+0.4fi ", s_potential_matrix[IDC2D(count1,count2,Nt)].x,s_potential_matrix[IDC2D(count1,count2,Nt)].y);
//    	}
//    	printf("\n");
//    }
//    error=cudaMemcpy(R,d_R,Nt*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
//    printf("the test upper triangular matrix in kernel is:\n");
//    for(count1=0;count1<Nt;count1++)
//       {
//       	for(int count2=0;count2<Nt;count2++)
//       	{
//       		printf("%0.4f%+0.4fi ", R[IDC2D(count1,count2,Nt)].x,R[IDC2D(count1,count2,Nt)].y);
//       	}
//       	printf("\n");
//       }

//    error=cudaMemcpy(Eu,d_Eu,pathNum*sizeof(float),cudaMemcpyDeviceToHost);
//    if(error!=cudaSuccess)
//    {
////   	printf("Eu returned error code %d, line %d\n", error, __LINE__);
////    			 	exit(EXIT_FAILURE);
//   	}


 //fine out the symbol vector index among all the block output Euclidean distance
//    int Eu_mini_index=0;
//    float Eu_mini_value=Eu[0];
//    for(count1=0;count1<pathNum;count1++)
//    {
//      if(Eu[count1]<Eu_mini_value)
//      {
//    	  Eu_mini_value=Eu[count1];
//    	  Eu_mini_index=count1;
//      }
//    }
    float *d_Eu;
    cudaMalloc((void**) &d_Eu,pathNum*sizeof(float));
    cudaMemcpy(d_Eu,Eu,pathNum*sizeof(float),cudaMemcpyHostToDevice);
    int *Eu_mini_index=(int*)malloc(sizeof(int));
    cublasIsamin(handle,pathNum,d_Eu,1,Eu_mini_index);
    if(*Eu_mini_index-1<pathNum/4)
    {
    for(count1=0;count1<Nt;count1++)
    {
     s_kernel[list[count1]-1]=s_potential1[IDC2D((MATRIX_SIZE-count1-1),(*Eu_mini_index-1),pathNum/4)];
    }
    }
    else if(*Eu_mini_index-1<pathNum/2)
    {
        for(count1=0;count1<Nt;count1++)
        {
         s_kernel[list[count1]-1]=s_potential2[IDC2D((MATRIX_SIZE-count1-1),(*Eu_mini_index-1-pathNum/4),pathNum/4)];
        }
    }
    else if(*Eu_mini_index-1<(pathNum*3)/4)
    {
        for(count1=0;count1<Nt;count1++)
        {
         s_kernel[list[count1]-1]=s_potential3[IDC2D((MATRIX_SIZE-count1-1),(*Eu_mini_index-1-pathNum/2),pathNum/4)];
        }
    }
    else
    {
        for(count1=0;count1<Nt;count1++)
        {
         s_kernel[list[count1]-1]=s_potential3[IDC2D((MATRIX_SIZE-count1-1),(*Eu_mini_index-1-(pathNum*3)/4),pathNum/4)];
        }
    }

//    	for(int count2=0;count2<MATRIX_SIZE;count2++)
//    	{
//    		printf("%0.4f%+0.4f ", s_share_matrix[IDC2D(count1,count2,MATRIX_SIZE)].x, s_share_matrix[IDC2D(count1,count2,MATRIX_SIZE)].y);
//    	}
//    	printf("\n");
//    }
//       error=cudaMemcpy(s_hat,d_s_hat,Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
          if(error!=cudaSuccess)
          			   	{
//          			  	printf("s_hat returned error code %d, line %d\n", error, __LINE__);
          //			 	exit(EXIT_FAILURE);
          			   	}
//          for(count1=0;count1<pathNum;count1++)
//          {
//          printf("the unconstrained estimation is:\n");
//          	for(int count2=0;count2<MATRIX_SIZE;count2++)
//          	{
//          		printf("%0.4f%+0.4fi ", s_hat[count2].x, s_hat[count2].y);
//          	}
//          	printf("\n");
//          }

				ret=cublasDestroy(handle);
			    if (ret != CUBLAS_STATUS_SUCCESS)
			    {
//			        printf("cublasDestroy returned error code %d, line(%d)\n", ret, __LINE__);
//			        exit(EXIT_FAILURE);
			    }
			    printf("the s_kernel is :\n");
			    for(count1=0;count1<Nt;count1++)
			    {
			    	printf("%0.4f%+0.4fi ", s_kernel[count1].x, s_kernel[count1].y);
			    }
			    printf("\n");

			   	free(s_sub_index);
//			   	cudaFree(d_s_sub_index);
			   	cudaFree(d_list);
			   	free(Eu);
			   	cudaFree(d_Eu);
//			   	free(s_potential_matrix);
//			   	cudaFree(d_s_potential_matrix);
			   	cudaFree(d_psymbolconstellation);
			   	free(Eu_mini_index);
//			   	cudaFree(d_R);

//			    cuComplex *s_potential1,*s_potential2,*s_potential3,*s_potential4;
//			    cuComplex *d_s_potential1,*d_s_potential2,*d_s_potential3,*d_s_potential4;
//			 	float *Eu1,*Eu2,*Eu3,*Eu4;
//			 	float *d_Eu1,*d_Eu2,*d_Eu3,*d_Eu4;
//			 	int *sub_index1,*sub_index2,*sub_index3,*sub_index4;
//			 	int *d_sub_index1,*d_sub_index2,*d_sub_index3,*d_sub_index4;
			 	cudaFreeHost(s_potential1);cudaFreeHost(s_potential2);cudaFreeHost(s_potential3);cudaFreeHost(s_potential4);
			 	cudaFree(d_s_potential1);cudaFree(d_s_potential2);cudaFree(d_s_potential3);cudaFree(d_s_potential4);
			 	cudaFreeHost(Eu1);cudaFreeHost(Eu2);cudaFreeHost(Eu3);cudaFreeHost(Eu4);
			 	cudaFree(d_Eu1);cudaFree(d_Eu2);cudaFree(d_Eu3);cudaFree(d_Eu4);
			 	cudaFreeHost(sub_index1);cudaFreeHost(sub_index2);cudaFreeHost(sub_index3);cudaFreeHost(sub_index4);
			 	cudaFree(d_sub_index1);cudaFree(d_sub_index2);cudaFree(d_sub_index3);cudaFree(d_sub_index4);
}

