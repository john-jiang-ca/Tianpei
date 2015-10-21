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
#define blockNum 4
#define stride 1
__global__ void FEpath(
		cuComplex *R,  //upper triangular matrix after cholesky factorization
		cuComplex *s_hat,  //unconstrained estimation of transmitted symbol vector s
//		cuComplex *s_matrix_share,
		cuComplex *s_potential_matrix,   //the matrix use to store all the solution candidates from all the blocks
		int *s_sub_index,   //full factorial index matrix
//		cuComplex *s,  //decoding results
		float *Eu,  //Euclidean distance
		int pitch_R,    //the number of transmit antennas
		int pitch_index,    //the number of receive antennas
		int pitch_p,
		int M,    //modulation scheme
		int threadNum,    //number of threads
		int *list,     //the permutation list
		cuComplex *psymbolconstellation,
		int index





		)
{
	//need to consider the resource allocation
	int tx=blockIdx.x*blockDim.x+threadIdx.x;     //if the path number is small we can allocate the kernel into one block so that we can use the shared memory
int tid=threadIdx.x;
int Nt=MATRIX_SIZE;
//int bid=blockIdx.x;
//allocate shared memory
	extern __shared__ cuComplex array[];
//	extern __shared__ cuComplex array2[];
//	  extern __shared__ float Eu_vector[];
//	 extern __shared__ float Euclidean[];
//    __shared__ int mini_Eu_index_temp;
//    __shared__ int mini_Eu_index;
//    extern __shared__ int s_sub_index[];
	error_t error;
	int count1, count2,count3,count4;
	__shared__ float d;    //the minimum distance unit between the signal constellation, the distance is usually 2d
	int rho=ceil(sqrt(float(Nt))-1);
	int pathNum=pow(float(M),float(rho));


//	int blockNum=pathNum/(threadNum);
//	int threadNum=1024;free(R_Eu_share);

//	float *distance;
//	int *resu
//#if __CUDA_ARCH__ >+300

//#endif
	   __shared__ cuComplex alpha, beta;
	   alpha.x=1;
	   alpha.y=0;
	   beta.x=0;
	   beta.y=0;

//    float *Eu_vector=Euclidean;
	cuComplex *R_share=array+Nt+2*threadNum;   //the upper trian
//	cuComplex *s_matrix_share=array+Nt*Nt; //the full expansion detection different for different path
//	cuComplex *s_sub_share=array+Nt*(Nt+1)/2;
//	cuComplex *s_sub_share=(cuComplex*)malloc(sizeof(cuComplex))
//    cuComplex *s_hat_share=array+Nt*(Nt+1)/2;
	cuComplex *s_hat_share=array;
    cuComplex *s_temp=array+Nt;
//    cuComplex s_temp[tid];
//    cuComplex *R_Eu_share=array+Nt*(Nt+1)/2+Nt;
    cuComplex *R_Eu_share=(cuComplex*)malloc(Nt*sizeof(cuComplex));
    cuComplex *Eu_norm_share=array+Nt+threadNum;
//    cuComplex Eu_norm_share;
//   float  *Eu_vector=(float*)malloc(threadNum*sizeof(float));
//    cuComplex *s_mini_share=array+Nr*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum+Nt*threadNum+threadNum;


//	   __syncthreads();
    //single expansion
    //
//    for(int i=0;i<pathNum;i+=threadNum)
//    {
//    if((threadNum*blockNum*index+tx)<pathNum)
//    {
   if(tid==0)
   {
    for(count1=0;count1<MATRIX_SIZE;count1++)
    {

		for(count2=0;count2<MATRIX_SIZE;count2++)
		{
//			R_share[(2*MATRIX_SIZE-tid+1)*(tid)/2+count2]=R[IDC2D((count2+tid),tid,pitch_R)];
			R_share[IDC2D(count2,count1,pitch_R)]=R[IDC2D(count2,count1,pitch_R)];
		}


		s_hat_share[count1]=s_hat[count1];

//		printf("the s_hat_share is unchanged? %d, %0.4f%+0.4fi\n", tid, s_hat_share[tid]);
    }
   }
	__syncthreads();
//
//	if(tx==0)
//	{
//		printf("the original upper triangular matrix before decoding is:\n");
//		  for(count1=0;count1<MATRIX_SIZE;count1++)
//		  {
//			  for(count2=0;count2<MATRIX_SIZE;count2++)
//			  {
//				  printf("%0.4f%+0.4fi ", R[IDC2D(count1,count2,pitch_R)].x,R[IDC2D(count1,count2,pitch_R)].y);
//
//			  }
//			  printf("\n");
//		  }
//
//			printf("the shared upper triangular matrix before decoding is:\n");
//			  for(count1=0;count1<MATRIX_SIZE;count1++)
//			  {
//				  for(count2=0;count2<MATRIX_SIZE;count2++)
//				  {
//					  printf("%0.4f%+0.4fi ", R_share[IDC2D(count1,count2,pitch_R)].x,R_share[IDC2D(count1,count2,pitch_R)].y);
//
//				  }
//				  printf("\n");
//			  }
////			  printf("the sub_index is\n");
//////				printf("the shared upper triangular matrix before decoding is:\n");
////				  for(count1=pathNum-16;count1<pathNum;count1++)
////				  {
////					  for(count2=0;count2<rho;count2++)
////					  {
////						  printf("%d ", s_sub_index[IDC2D(count2,count1,pitch_index)]);
////
////					  }
////					  printf("\n");
////				  }
//	}
//		  printf("s_hat is!!!!\n");
//		  for(count1=0;count1<MATRIX_SIZE;count1++)
//		  {
//			  printf("%0.4f%+0.4fi ",s_hat_share[count1].x,s_hat_share[count1].y );
//		  }
//		  printf("kernel is working\n");
//	}
//	for(count4=0;count4<5;count4++)
//	{
//	for(count1=0;count1<rho;count1++)
//	{
//	s_sub_share[IDC2D(tid,count1,rho)]=psymbolconstellation[s_sub_index[IDC2D((index*blockNum*threadNum+tx),count1,rho)]];//read by column major so that it can be continue
//	}
//	if(tx==0)
//	{
//	printf("s_hat_share is\n");
////	for(count2=0;count2<blockDim.x;count2++)
////	{
//	for(count1=0;count1<MATRIX_SIZE;count1++)
//	{
//	printf("%0.4f%+0.4fi ", s_hat_share[count1].x,s_hat_share[count1].y);
//	}
//	printf("\n");
////	}
//	}
//Eu_vector[tid]=0;
	Eu[blockNum*threadNum*index+tx]=0;
for (count1=Nt-1; count1>=0; count1--)
{

		if (count1<Nt-rho)
		{
			s_temp[tid]=s_hat_share[count1];
//			if(tx==0)
//			{
//				printf("the s_temp[tid] is not wrong!!%0.4f%+0.4fi: ", s_hat_share[count1].x,s_hat_share[count1].y);
//				printf("\n");
//			}
			for (count2=count1+1;count2<Nt; count2++)
			{
				s_temp[tid]=complex_add(s_temp[tid],complex_mulcom(complex_div(R_share[IDC2D(count2,count1,pitch_R)],R_share[IDC2D(count1,count1,pitch_R)]),(complex_sub(s_hat_share[count2],s_potential_matrix[IDC2D(count2,(index*blockNum*threadNum+tx),pitch_p)]))));
			}
//			if(tx==0)
//			{
//
//				printf("now the s_temp[tid] is %0.4f%+0.4fi: ", s_temp[tid].x,s_temp[tid].y);
//			}
//			s_share[IDC2D(tid,count1,Nt)]=R_temp[tid];

	//mapping the estimation s_share to the correlated signal constellation(to be continued)
	if(M==2)   //BPSK
	{
//     for(count1=0;count1<Nt-rho; count1++)
//     {
    	 d=sqrt(float(float(1)/float(Nt)));
//    	 s_share_Q[IDC2D(threadIdx.x,count1,Nt)].x=float(int(R_temp[threadIdx.x].x/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//    	 s_share_Q[IDC2D(threadIdx.x,count1,Nt)].y=0;
//     }0.6565-1.0724i 0.1581+0.1581i 0.3595-0.7887i 0.1581-0.4743i

    	 if(s_temp[tid].x>0)
    	 {
    		 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=d;
    		 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=0;
    	 }
    	 else
    	 {
    		 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=(-d);
    		 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=0;
    	 }
	}
	else if(M==4)   //QPSK
	{
//		if(tx==0)//		printf("the up triangular matrix is:\n");
		//		  for(count1=0;count1<MATRIX_SIZE;count1++)
		//		  {
		//			  for(count2=0;count2<MATRIX_SIZE;count2++)
		//			  {
		//				  printf("%0.4f%+0.4fi ", R_share[IDC2D(count1,count2,MATRIX_SIZE)].x,R_share[IDC2D(count1,count2,MATRIX_SIZE)].y);
		//
		//			  }
		//			  printf("\n");
		//		  }
		//		  printf("kernel is working\n");
//		{
//       #if __CUDA_ARCH__ >=300
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
        	   distance[count2]=sqrt(powf(s_temp[tid].x-(-d),2)+pow(s_temp[tid].y-0,2)); break;
        	   case 1:
        	   distance[count2]=sqrt(powf(s_temp[tid].x-0,2)+pow(s_temp[tid].y-(-d),2)); break;
        	   case 2:
        	   distance[count2]=sqrt(pow(s_temp[tid].x-(d),2)+pow(s_temp[tid].y-0,2)); break;
        	   case 3:
        	   distance[count2]=sqrt(pow(s_temp[tid].x-0,2)+pow(s_temp[tid].y-d,2)); break;
               default:
                   #if __CUDA_ARCH__ >=300
            	   printf("result error code %d\n", error);
                   #endif
            	   break;
        	   }
           }
//#if __CUDA_ARCH__ >=300
//           ret=cublasIsamin(handle,4,distance,1,result);
//#endif
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



//           {
//             #if __CUDA_ARCH__ >=300
//        	   printf("cublasIsamin failed return error code %d, line %d\n",error,__LINE__);
////        	   exit(EXIT_bbbb
//#endif
//           }

           switch (mini_index)
           {
           case 1:
        	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=-d;
        	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=0;
        	break;
           case 2:
           	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=0;
           	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=-d;
           	break;
           case 3:
           	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=d;
           	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=0;
           	break;
           case 4:
       	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=0;
       	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=d;
       	  break;
           default:
               #if __CUDA_ARCH__ >=300
        	   printf("result error code %d\n", error);
               #endif
           }
           free(distance);
           free(result);
	}
//#if __CUDA_ARCH__ >=300
//          cublasDestroy(handle);
//#endif
	else if(M==16)  //16QAM
	{
//	     for(count1=0;count1<Nt-rho; count1++)
//	     {
	    	 d=sqrt(float(3)/(2* (float)(Nt*(M-1))));
//	    	 if(tx==0)
//	    	 {
//	    	 printf("the distance of constellation is: %0.4f ", d);
//	    	 printf("the s_temp[tid] is %0.4f%+0.4fi ", s_temp[tid].x,s_temp[tid].y);
//	    	 }
//	  s_share_Q[IDC2D(tid,count1,Nt)].x=float(int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//	  s_share_Q[IDC2D(threadIdx.x,count1,Nt)].y=float(int(s_share[IDC2D(threadIdx.x,count1,Nt)].y/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//	     if(s_temp[tid].x<(-2*d))
//	    	{
//	    	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=(-3*d);
//	    	}
//	    else if(s_temp[tid].x>(2*d))
//	    	{
//	    	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=(3*d);
//	    	}
//	    else if(s_temp[tid].x>=0&&s_temp[tid].x<=2*d)
//	    	{
//	    	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=d;
//	    	}
//	    else if(s_temp[tid].x>=(-2*d)&&s_temp[tid].x<=0)
//	        {
//	    	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=(-d);
//	    	}
//
//
//
//
//	    if(s_temp[tid].y<(-2*d))
//	    	{
//	    	   s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=(-3*d);
//	    	}
//	   else if(s_temp[tid].y>(2*d))
//	    	 {
//	    	  s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=(3*d);
//	         }
//	   else if(s_temp[tid].y>=0&&s_temp[tid].y<=(2*d))
//	        {
//	    	s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=d;
//	    	}
//	    else if(s_temp[tid].y>=(-2*d)&&s_temp[tid].y<=0)
//	       {
//	        s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=(-d);
//	       }
//	     }

	    	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].x=2*d*int(s_temp[tid].x/(2*d))+d*(int(s_temp[tid].y/sqrt(pow(s_temp[tid].x,2)+pow(s_temp[tid].x,2))));
	    	 s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)].y=2*d*int(s_temp[tid].y/(2*d))+d*(int(s_temp[tid].y/sqrt(pow(s_temp[tid].y,2)+pow(s_temp[tid].y,2))));
	}
	else if(M==64)   //64QAM
		{
//	     for(count1=0;count1<Nt-rho; count1++)
//	     {
	    	 d=sqrt(3/(2* (float)(Nt*(M-1))));
//	     s_potential_matrix[IDC2D(tid,count1,Nt)].x=float(int(s_potential_matrix[IDC2D(threadIdx.x,count1,Nt)].x/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//	    s_matrix_Q[IDC2D(threadIdx.x,count1,Nt)].y=float(int(s_potential_matrix[IDC2D(threadIdx.x,count1,Nt)].y/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
	     }
//		}
		}
	else
	{
//			s_share[IDC2D(tid,count1,Nt)]=s_sub_share[IDC2D(tx,(Nt-count1-1),Nt)];
//		s_potential_matrix[IDC2D(count2,(index*blockNum*threadNum+tx),pitch_p)]=s_sub_share[IDC2D(tid,(Nt-count1-1),rho)];
		s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)]=psymbolconstellation[s_sub_index[IDC2D((Nt-count1-1),(index*blockNum*threadNum+tx),pitch_index)]];
//		printf("the s_potential_matrix in each thread is:\n");
//		printf("%0.4f%+0.4fi ", s_potential_matrix[IDC2D(tid,count1,Nt)].x,s_potential_matrix[IDC2D(tid,count1,Nt)].y);
//		printf("\n");
	}
//		if(tid==0)
//		{
//
//			printf("the decoding matrix is:%0.4f%+0.4fi ",s_potential_matrix[IDC2D(count2,(index*blockNum*threadNum+tx),pitch_p)] );
//		}

//		R_Eu_share[IDC2D(tid,count1,Nt)]=complex_sub(s_potential_matrix[IDC2D(count2,(index*blockNum*threadNum+tx),pitch_p)],s_hat_share[count1]);
		R_Eu_share[count1]=complex_sub(s_potential_matrix[IDC2D(count1,(index*blockNum*threadNum+tx),pitch_p)],s_hat_share[count1]);
		Eu_norm_share[tid]=beta;
		for(count3=count1;count3<MATRIX_SIZE;count3++)
		{
		Eu_norm_share[tid]=complex_add(Eu_norm_share[tid],complex_mulcom(R_share[IDC2D(count3,count1,pitch_R)],R_Eu_share[count3]));
		}
//		Eu_vector[tid]=Eu_vector[tid]+pow(Eu_norm_share.x,2)+pow(Eu_norm_share.y,2);
		Eu[(index*blockNum*threadNum+tx)]=Eu[(index*blockNum*threadNum+tx)]+pow(Eu_norm_share[tid].x,2)+pow(Eu_norm_share[tid].y,2);

	}
free(R_Eu_share);
//free(Eu_vector);
__syncthreads();
//if(tid>=0&&tid<threadNum)
//{
//	Eu[(index*blockNum*threadNum+tx)]=Eu_vector[tid];
//}
//__syncthreads();
//    }

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
//	cuComplex *ss;
//	ss=(cuComplex*)malloc(MATRIX_SIZE*sizeof(cuComplex));
//	cuComplex *s_sub;
//	s_sub=(cuComplex*)malloc(pow(M,rho)*rho*sizeof(cuComplex));   //all the possible full expansion sub vector
	int  pathNum;
	pathNum=pow(M,rho);
	int *d_list,*d_s_sub_index;
	int *s_sub_index=(int*)calloc(1,rho*pow(M,rho)*sizeof(int));
	fullfact(rho,M,s_sub_index);    //get  the indexes of all the possible rho length symbol vectors
//	int blockNum=BLOCK_NUM;   //determined by the path number
//	int pathNum=pow(M,rho);  //number of search path
	int threadNum=ceil(pathNum/(blockNum*stride)); //determined by the path number
	float *Eu,*d_Eu;
	Eu=(float*)calloc(1,blockNum*sizeof(float));
	cuComplex *s_potential_matrix=(cuComplex*)calloc(1,pathNum*Nt*sizeof(cuComplex));
	cuComplex  *d_R, *d_s_potential_matrix,*d_psymbolconstellation;
	int count1;
//	int *j;
//	j=(int*)malloc(sizeof(int));
	cublasHandle_t handle;
		cublasStatus_t ret;
		cudaError_t error;
		size_t pitch_R,pitch_potential,pitch_index;
		ret=cublasCreate(&handle);
	    error=cudaMallocPitch((void**) &d_R, &pitch_R, MATRIX_SIZE*sizeof(cuComplex),MATRIX_SIZE);
	    error=cudaMallocPitch((void**) &d_s_sub_index, &pitch_index,pathNum*sizeof(int),rho);
		error=cudaMallocPitch((void**) &d_s_potential_matrix,&pitch_potential, pathNum*sizeof(cuComplex),Nt);
		error=cudaMalloc((void**) &d_list, MATRIX_SIZE*sizeof(int));
		error=cudaMalloc((void**) &d_Eu, pathNum*sizeof(float));
        error=cudaMalloc((void**) &d_psymbolconstellation, M*sizeof(cuComplex));
        error=cudaMemcpy2D(d_R,pitch_R,R,sizeof(cuComplex),sizeof(cuComplex),MATRIX_SIZE*MATRIX_SIZE,cudaMemcpyDeviceToDevice);
		error=cudaMemcpy(d_psymbolconstellation, psymbolconstellation, M*sizeof(cuComplex),cudaMemcpyHostToDevice);
		error=cudaMemcpy2D(d_s_sub_index,pitch_index, s_sub_index,pathNum*sizeof(int), int(pow(M,rho))*sizeof(int),rho,cudaMemcpyHostToDevice);
		error=cudaMemcpy(d_list, list, Nt*(sizeof(int)),cudaMemcpyHostToDevice);
	 int sharedMem;
    sharedMem=6000*sizeof(cuComplex);
   float duration;
   clock_t start, end;
   start=clock();
   for(count1=0;count1<stride;count1++)
   {
	FEpath<<<blockNum, threadNum,sharedMem>>>(d_R, d_s_hat, d_s_potential_matrix,d_s_sub_index, d_Eu, int(pitch_R/sizeof(cuComplex)), int(pitch_index/sizeof(int)),int (pitch_potential/sizeof(cuComplex)), M,threadNum,d_list,d_psymbolconstellation,count1);
	error=cudaDeviceSynchronize();
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
   }
	end=clock();
//	cudaProfilerStop();
	duration=double(end-start);
	printf("hey %0.4f ", duration);
	printf("\n");

		if(error!=cudaSuccess)
		{
//		printf("cudaDeviceSynchronize returned error code %d, line %d\n", error, __LINE__);
//				 	exit(EXIT_FAILURE);
		}

//	printf("Eu_num is %d", Eu_num);
//    error=cudaMemcpy(s_potential_matrix,d_s_potential_matrix,pathNum*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
error=cudaMemcpy2D(s_potential_matrix,pathNum*sizeof(cuComplex), d_s_potential_matrix,pitch_potential, pathNum*sizeof(cuComplex),MATRIX_SIZE,cudaMemcpyDeviceToHost);
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
    int *Eu_mini_index=(int*)malloc(sizeof(int));
    cublasIsamin(handle,pathNum,d_Eu,1,Eu_mini_index);
    for(count1=0;count1<Nt;count1++)
    {
     s_kernel[list[count1]-1]=s_potential_matrix[IDC2D((MATRIX_SIZE-count1-1),(*Eu_mini_index-1),Nt)];
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
			   	cudaFree(d_s_sub_index);
			   	cudaFree(d_list);
			   	free(Eu);
			   	cudaFree(d_Eu);
			   	free(s_potential_matrix);
			   	cudaFree(d_s_potential_matrix);
			   	cudaFree(d_psymbolconstellation);
			   	free(Eu_mini_index);
//			   	cudaFree(d_R);

}

