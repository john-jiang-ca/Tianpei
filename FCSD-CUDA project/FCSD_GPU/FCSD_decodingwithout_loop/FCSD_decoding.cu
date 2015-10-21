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
#define BLOCK_NUM 64
#define stride 256
__global__ void FEpath(
		cuComplex *R,  //upper triangular matrix after cholesky factorization
		cuComplex *s_hat,  //unconstrained estimation of transmitted symbol vector s
		cuComplex *s_potential_matrix,   //the matrix use to store all the solution candidates from all the blocks
		int *s_sub_index,   //full factorial index matrix
//		cuComplex *s,  //decoding results
		float *Eu,  //Euclidean distance
		int Nt,    //the number of transmit antennas
		int Nr,    //the number of receive antennas
		int M,    //modulation scheme
		int threadNum,    //number of threads
		int *list,     //the permutation list
		cuComplex *psymbolconstellation




		//the variable for test
		)//    cuComplex *Eu_norm_mini_share=array+Nr*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum+Nt*threadNum+threadNum;
{
	//need to consider the resource allocation
	int tx=blockIdx.x*blockDim.x+threadIdx.x;     //if the path number is small we can allocate the kernel into one block so that we can use the shared memory
int tid=threadIdx.x;
int bid=blockIdx.x;
//allocate shared memory
	extern __shared__ cuComplex array[];
	 extern __shared__ float Eu_vector[];
    __shared__ int mini_Eu_index_temp;
    __shared__ int mini_Eu_index;
//    extern __shared__ int s_sub_index[];
	error_t error;
	int count1, count2,count3,count4,lag;
	__shared__ float d;    //the minimum distance unit between the signal constellation, the distance is usually 2d
	int rho=ceil(sqrt(float(Nt))-1);
	int pathNum=pow(float(M),float(rho));


//	int blockNum=pathNum/(threadNum);
//	int threadNum=1024;

//	float *distance;
//	int *resu
//#if __CUDA_ARCH__ >+300

//#endif
	   __shared__ cuComplex alpha, beta;
	   alpha.x=1;
	   alpha.y=0;
	   beta.x=0;
	   beta.y=0;
//	cuComplex *R_Eu_share;
//	   if(tx==1)
//	   {
/*R_share=Nr*Nt;
 * s_matrix_share: Nt*threadNumerror=cudaMemcpy(d_R, R, Nr*Nt*sizeof(cuComplex),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMemcpy d_R returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
 * s_sub_share: rho*threadNum
 * s_hat_share : Nt*1
 * s_temp: threadNum*1
 * R_Eu_share:Nt*threadNum
 * Eu_norm_share : threadNum*1
// * Eu_norm_mini_share :1
 * s_mini_share:  Nt*1
 */

	cuComplex *R_share=array;   //the upper trian
	cuComplex *s_matrix_share=array+Nt*Nt; //the full expansion detection different for different path
	cuComplex *s_sub_share=array+Nt*Nt+Nt*threadNum;
    cuComplex *s_hat_share=array+Nt*Nt+Nt*threadNum+rho*threadNum;
    cuComplex *s_temp=array+Nt*Nt+Nt*threadNum+rho*threadNum+Nt;
    cuComplex *R_Eu_share=array+Nt*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum;
    cuComplex *Eu_norm_share=array+Nt*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum+Nt*threadNum;
//    cuComplex *s_mini_share=array+Nr*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum+Nt*threadNum+threadNum;

	//in each thread there are one column sthore s_hat, s_share s_hat_share s_share_Q and s_sub_share
//	cuComplex *s_share=array+Nr*Nt+threadNum*rho+Nt;  //the single expansion before quantization
//	cuComplex *R_temp=array+Nr*Nt+threadNum*rho*+Nt+threadNum*Nt;
//	cuComplex *s_share_Q=array+Nr*Nt+threadNum*rho+(Nt+threadNum)+threadNum*Nt;//the single expansion after quantization
//    cuComplex *R_Eu_share1=array+Nr*Nt+thread
 //used to store the vector of R(s_Q-s_hat);
//    cuComplex *R_Eu_share2=array+Nr*Nt+threadNum*rho+(Nt+threadNum)+threadNum*Nt+2*threadNum*Nt;//stored in column major
//    float *distance=array+Nr*Nt*sizeof(cuComplex)+pathNum*rho*sizeof(cuComplex)+(Nt+1)*sizeof(cuComplex)+2*pathNum*Nt*sizeof(cuComplex)+Nr*pathNum*sizeof(cuComplex);
//	   }
//	   __syncthreads();
    //single expansion
    //
//    for(int i=0;i<pathNum;i+=threadNum)
//    {
    if(tid>=0&&tid<MATRIX_SIZE)
    {
//	for(count1=0;count1<Nr;count1++)
//	{
		for(count2=0;count2<MATRIX_SIZE;count2++)
		{
			R_share[IDC2D(tid,count2,Nt)]=R[IDC2D(tid,count2,Nt)];
		}
//	}

		s_hat_share[tid]=s_hat[tid];
    }
	__syncthreads();
//	if(tx==0)
//	{
//		printf("the upper triangular matrix before decoding is:\n");
//		  for(count1=0;count1<MATRIX_SIZE;count1++)
//		  {
//			  for(count2=0;count2<MATRIX_SIZE;count2++)
//			  {
//				  printf("%0.4f%+0.4fi ", R_share[IDC2D(count1,count2,MATRIX_SIZE)].x,R_share[IDC2D(count1,count2,MATRIX_SIZE)].y);
//
//			  }
//			  printf("\n");
//		  }
//		  printf("s_hat is!!!!\n");
//		  for(count1=0;count1<MATRIX_SIZE;count1++)
//		  {
//			  printf("%0.4f%+0.4fi ",s_hat_share[count1].x,s_hat_share[count1].y );
//		  }
//		  printf("kernel is working\n");
//	}
	for(count4=0;count4<5;count4++)
	{
	for(count1=0;count1<rho;count1++)
	{
	s_sub_share[IDC2D(tid,count1,rho)]=psymbolconstellation[s_sub_index[IDC2D(((blockIdx.x*blockDim.x*stride)+blockDim.x*count4+threadIdx.x),count1,rho)]];//read by column major so that it can be continue
	}
for (count1=Nt-1; count1>=0; count1--)
{
		if (count1<Nt-rho)
		{
			s_temp[tid]=s_hat_share[count1];
			for (count2=count1+1;count2<Nt; count2++)
			{
				s_temp[tid]=complex_add(s_temp[tid],complex_mulcom(complex_div(R_share[IDC2D(count1,count2,Nt)],R_share[IDC2D(count1,count1,Nt)]),(complex_sub(s_hat_share[count2],s_matrix_share[IDC2D(tid,count2,Nt)]))));
			}
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
    		 s_matrix_share[IDC2D(tid,count1,Nt)].x=d;
    		 s_matrix_share[IDC2D(tid,count1,Nt)].y=0;
    	 }
    	 else
    	 {
    		 s_matrix_share[IDC2D(tid,count1,Nt)].x=(-d);
    		 s_matrix_share[IDC2D(tid,count1,Nt)].y=0;
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
        	   distance[count2]=sqrt(pow(s_temp[tid].x-(-d),2)+pow(s_temp[tid].y-0,2)); break;
        	   case 1:
        	   distance[count2]=sqrt(pow(s_temp[tid].x-0,2)+pow(s_temp[tid].y-(-d),2)); break;
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
        	s_matrix_share[IDC2D(tid,count1,Nt)].x=-d;
        	s_matrix_share[IDC2D(tid,count1,Nt)].y=0;
        	break;
           case 2:
           	s_matrix_share[IDC2D(tid,count1,Nt)].x=0;
           	s_matrix_share[IDC2D(tid,count1,Nt)].y=-d;
           	break;
           case 3:
           	s_matrix_share[IDC2D(tid,count1,Nt)].x=d;
           	s_matrix_share[IDC2D(tid,count1,Nt)].y=0;
           	break;
           case 4:
       	  s_matrix_share[IDC2D(tid,count1,Nt)].x=0;
       	  s_matrix_share[IDC2D(tid,count1,Nt)].y=d;
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
	    	 d=sqrt(3/(2* (float)(Nt*(M-1))));
//	  s_share_Q[IDC2D(tid,count1,Nt)].x=float(int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//	  s_share_Q[IDC2D(threadIdx.x,count1,Nt)].y=float(int(s_share[IDC2D(threadIdx.x,count1,Nt)].y/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
	     if(s_temp[tid].x<(-2*d))
	    	{
	    	  s_matrix_share[IDC2D(tid,count1,Nt)].x=(-3*d);
	    	}
	    else if(s_temp[tid].x>(2*d))
	    	{
	    	  s_matrix_share[IDC2D(tid,count1,Nt)].x=(3*d);
	    	}
	    else if(s_temp[tid].x>=0&&s_temp[tid].x<=2*d)
	    	{
	    	  s_matrix_share[IDC2D(tid,count1,Nt)].x=d;
	    	}
	    else if(s_temp[tid].x>=(-2*d)&&s_temp[tid].x<=0)
	        {
	    	  s_matrix_share[IDC2D(tid,count1,Nt)].x=(-d);
	    	}




	    if(s_temp[tid].y<(-2*d))
	    	{
	    	   s_matrix_share[IDC2D(tid,count1,Nt)].y=(-3*d);
	    	}
	   else if(s_temp[tid].y>(2*d))
	    	 {
	    	  s_matrix_share[IDC2D(tid,count1,Nt)].y=(3*d);
	         }
	   else if(s_temp[tid].y>=0&&s_temp[tid].y<=(2*d))
	        {
	    	s_matrix_share[IDC2D(tid,count1,Nt)].y=d;
	    	}
	    else if(s_temp[tid].y>=(-2*d)&&s_temp[tid].y<=0)
	       {
	        s_matrix_share[IDC2D(tid,count1,Nt)].y=(-d);
	       }
//	     }
	}
	else if(M==64)   //64QAM
		{
//	     for(count1=0;count1<Nt-rho; count1++)
//	     {
	    	 d=sqrt(3/(2* (float)(Nt*(M-1))));
//	     s_matrix_share[IDC2D(tid,count1,Nt)].x=float(int(s_matrix_share[IDC2D(threadIdx.x,count1,Nt)].x/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
//	    s_matrix_Q[IDC2D(threadIdx.x,count1,Nt)].y=float(int(s_matrix_share[IDC2D(threadIdx.x,count1,Nt)].y/d)*d)+float(int(round(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d-int(s_share[IDC2D(threadIdx.x,count1,Nt)].x/d)))*d);
	     }
//		}
		}
	else
	{
//			s_share[IDC2D(tid,count1,Nt)]=s_sub_share[IDC2D(tx,(Nt-count1-1),Nt)];
		s_matrix_share[IDC2D(tid,count1,Nt)]=s_sub_share[IDC2D(tid,(Nt-count1-1),rho)];
//		printf("the s_matrix_share in each thread is:\n");
//		printf("%0.4f%+0.4fi ", s_matrix_share[IDC2D(tid,count1,Nt)].x,s_matrix_share[IDC2D(tid,count1,Nt)].y);
//		printf("\n");
	}
	}
//__syncthreads();


//calculation of the Euclidean distance of all the symbol vectors in one block
Eu_vector[tid]=0;
for(count2=0;count2<Nt;count2++)
{
R_Eu_share[IDC2D(tid,count2,Nt)]=complex_sub(s_matrix_share[IDC2D(tid,count2,Nt)],s_hat_share[count2]);
}

for(count2=0;count2<Nt;count2++)
{
Eu_norm_share[tid]=beta;
for(count3=0;count3<Nt;count3++)
{
Eu_norm_share[tid]=complex_add(Eu_norm_share[tid],complex_mulcom(R_share[IDC2D(count2,count3,Nt)],R_Eu_share[IDC2D(tid,count2,Nt)]));
}
		Eu_vector[tid]=Eu_vector[tid]+pow(Eu_norm_share[tid].x,2)+pow(Eu_norm_share[tid].y,2);
}

__syncthreads();
//find out the symbol vector candidate with minimum Euclidean distance in one block
if(tid==0)
{
//	printf("let go Eu_vector!!\n");
//	for(count1=0;count1<MATRIX_SIZE;count1++)
//	{
//		printf("%0.4f ", Eu_vector[count1]);
//	}
       __shared__ float  mini_Eu;
	     mini_Eu=Eu_vector[0];
         mini_Eu_index=0;
	     for(count3=0;count3<threadNum;count3++)
	     {
	    	if(Eu_vector[count3]<mini_Eu)
	    	{
	    		mini_Eu=Eu_vector[count3];
	    		mini_Eu_index_temp=count3;
	    	}
	     }
lag=0;
if(count4==0)
{
Eu[bid]=mini_Eu;
mini_Eu_index=mini_Eu_index_temp;
lag=1;
}
else
{
if(mini_Eu<Eu[bid])
{
	Eu[bid]=mini_Eu;
	mini_Eu_index=mini_Eu_index_temp;
	lag=1;
}
}

}
//__syncthreads();
if(tid>=0&&tid<=(MATRIX_SIZE-1))
{
if(lag==1)
{
	s_potential_matrix[IDC2D(bid,(list[tid]-1),Nt)]=s_matrix_share[IDC2D(mini_Eu_index,(MATRIX_SIZE-tid-1),Nt)];
}
//	s_potential_matrix[IDC2D(bid,(list[tid]-1),Nt)].x=1;
//	s_potential_matrix[IDC2D(bid,(list[tid]-1),Nt)].y=1;
}
__syncthreads();

	}
//}
//if(tid==0)
//{
////	printf("the s_potential_matrix is: \n");
////	for(count1=0;count1<threadNu;count1++)	int stride=512;
////	{
//
////printf("the s_share_matrix after iteration is:\n");
////for(count1=0;count1<threadNum;count1++)
////{
////	for(count2=0;count2<Nt;count2++)
////	{
////	printf("%0.4f%+0.4fi ", s_matrix_share[IDC2D(count1,count2,Nt)].x, s_matrix_share[IDC2D(count1,count2,Nt)].y);
////	}
////	printf("\n");
////}
//
//		for(count2=0;count2<Nt;count2++)
//		{
//		printf("the %d th row of the s_potential_matrix is :%0.4f%+0.4fi ", bid, s_potential_matrix[IDC2D(bid,count2,Nt)].x, s_potential_matrix[IDC2D(bid,count2,Nt)].y);
//		}
//		printf("\n");
////	}
//	printf("the miniEuindex is:\n");//error=cudaMemcpy(pR, d_pR, Nt*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
//	//      	    if(error!=cudaSuccess)
//	//      	     {
//	//      	    	 printf("cudaMemcpy pR returned error code %d, line %d\n", error, __LINE__);
//	//      	    	  exit(EXIT_FAILURE);
//	//      	      }
//	printf("%d ", mini_Eu_index);
//	printf("the Eu_vector is:\n");
//	for(count1=0;count1<threadNum;count1++)
//	{
//		printf("%0.4f ", Eu_vector[count1]);
//	}
//	printf("\n");
//
//			printf("the upper triangular matrix after decoding is:\n");
//			  for(count1=0;count1<MATRIX_SIZE;count1++)
//			  {
//				  for(count2=0;count2<MATRIX_SIZE;count2++)
//				  {
//					  printf("%0.4f%+0.4fi ", R_share[IDC2D(count1,count2,MATRIX_SIZE)].x,R_share[IDC2D(count1,count2,MATRIX_SIZE)].y);
//
//				  }
//				  printf("\n");
//			  }
//			  printf("FEpath kernel is working\n");
//}

}
//	printf("the FEpath is working\n");
//	for(count1=0;count1<MATRIX_SIZE;count1++)
//	{
//		for(count2=0;count2<MATRIX_SIZE;count2++)
//		{
//			R[IDC2D(count1,count2,Nt)].x=1;
//			R[IDC2D(count1,count2,Nt)].y=1;
//		}
//	}//error=cudaMemcpy(pR, d_pR, Nt*Nt*sizeof(cuComplex),cudaMemcpyDeviceToH//	if(tid==0)
//	{ost);
//      	    if(error!=cudaSuccess)
//      	     {
//      	    	 printf("cudaMemcpy pR returned error code %d, line %d\n", error, __LINE__);
//      	    	  exit(EXIT_FAILURE);
//      	      }
//}
//__syncthreads();



//if(tx==0)
//{
//printf("the detected symbol vector is:\n");
//for(count1=0;count1<MATRIX_SIZE;count1++)
//{
//	printf("%0.4f%+0.4f ", s[count1].x, s[count1].y);
//}//error=cudaMemcpy(pR, d_pR, Nt*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
//      	    if(error!=cudaSuccess)
//      	     {
//      	    	 printf("cudaMemcpy pR returned error code %d, line %d\n", error, __LINE__);
//      	    	  exit(EXIT_FAILURE);
//      	      }
//}
//}

//host
void FCSD_decoding(
		cuComplex *d_R,  //upper triangular matrix after cholesky factorization store in device side
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
	int blockNum=BLOCK_NUM;   //determined by the path number
//	int pathNum=pow(M,rho);  //number of search path
	int threadNum=pathNum/(blockNum*stride); //determined by the path number
	float *Eu,*d_Eu;
	Eu=(float*)calloc(1,blockNum*sizeof(float));
	cuComplex *s_potential_matrix=(cuComplex*)calloc(1,blockNum*Nt*sizeof(cuComplex));
	cuComplex  *d_s_potential_matrix,*d_psymbolconstellation;
	int count1;
//	int *j;
//	j=(int*)malloc(sizeof(int));
	cublasHandle_t handle;
		cublasStatus_t ret;
		cudaError_t error;
		ret=cublasCreate(&handle);
	    if (ret != CUBLAS_STATUS_SUCCESS)
	    {
//	        printf("cublasCreate returned error code %d, line(%d)\n", ret, __LINE__);
//	        exit(EXIT_FAILURE);
	    }
	    error=cudaMalloc((void**) &d_s_sub_index, rho*pow(M,rho)*sizeof(int));
	    if(error!=cudaSuccess)
	    {
//	    	printf("cudaMalloc d_s_sub_index returned error code %d, line %d\n", error, __LINE__);
//	    	exit(EXIT_FAILURE);
	    }
//	error=cudaMalloc((void**) &d_R, Nr*Nt*sizeof(cuComplex));
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMalloc d_R returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
//			error=cudaMalloc((void**) &d_s_sub, rho*sizeof(cuComplex));
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMalloc d_s_sub returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
//			error=cudaMalloc((void**) &d_s_hat, Nt*sizeof(cuComplex));
//			if(error!=cudaSuccess)
//			{
////				printf("cudaMalloc d_s_hat returned error code %d, line %d\n", error, __LINE__);
////				exit(EXIT_FAILURE);
//			}
			error=cudaMalloc((void**) &d_s_potential_matrix, blockNum*Nt*sizeof(cuComplex));
			if(error!=cudaSuccess)
			{
//				printf("cudaMalloc d_R returned error code %d,	 line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
			}
			error=cudaMalloc((void**) &d_list, MATRIX_SIZE*sizeof(int));
			if(error!=cudaSuccess)
			{
//				printf("cudaMalloc d_list returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
			}
			error=cudaMalloc((void**) &d_Eu, blockNum*sizeof(float));
			if(error!=cudaSuccess)
			{
//				printf("cudaMalloc d_Eu returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
			}
         error=cudaMalloc((void**) &d_psymbolconstellation, M*sizeof(cuComplex));
	       if(error!=cudaSuccess)
	          {
//	        	printf("cudaMalloc d_psymbolconstellation returned error code %d, line %d\n", error, __LINE__);
//		        exit(EXIT_FAILURE);
	           }

//			error=cudaMemcpy(d_R, R, Nr*Nt*sizeof(cuComplex),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMemcpy d_R returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
//			error=cudaMemcpy(d_s_sub, s_sub, rho*sizeof(cuComplex),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMemcpy d_s_sub returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
//			error=cudaMemcpy(d_s_hat, s_hat, Nt*sizeof(cuComplex),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
//				printf("cudaMemcpy d_s_hat returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
//			}
//			error=cudaMemcpy(d_s_potential_matrix, s_potential_matrix, blockNum*Nt*sizeof(cuComplex),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
////				printf("cudaMemcpy d_s_potential_matrix returned error code %d, line %d\n", error, __LINE__);
////				exit(EXIT_FAILURE);
//			}
//			error=cudaMemcpy(d_Eu, Eu, (blockNum)*sizeof(float),cudaMemcpyHostToDevice);
//			if(error!=cudaSuccess)
//			{
////				printf("cudaMemcpy d_Eu returned error code %d, line %d\n", error, __LINE__);
////				exi256t(EXIT_FAILURE);
//			}
			error=cudaMemcpy(d_psymbolconstellation, psymbolconstellation, M*sizeof(cuComplex),cudaMemcpyHostToDevice);
						if(error!=cudaSuccess)
						{
//							printf("cudaMemcpy d_psynbolconstellation returned error code %d, line %d\n", error, __LINE__);
//							exit(EXIT_FAILURE);
						}
		error=cudaMemcpy(d_s_sub_index, s_sub_index, rho*int(pow(M,rho))*sizeof(int),cudaMemcpyHostToDevice);
				if(error!=cudaSuccess)
					{
//						printf("cudaMemcpy d_psynbolconstellation returned error code %d, line %d\n", error, __LINE__);
//							exit(EXIT_FAILURE);
				}
				error=cudaMemcpy(d_list, list, Nt*(sizeof(int)),cudaMemcpyHostToDevice);
					if(error!=cudaSuccess)
					{
//					printf("cudaMemcpy d_s_kernel returned error code %d, line %d\n", error, __LINE__);
//					exit(EXIT_FAILURE);
				    }
			//add the psymbolconstellation
//				size_t heapsize;
//				heapsize=1024*sizeof(float);
				int sharedMem;
//				sharedMem=(Nr*Nt+Nt*threadNum+rho*threadNum+Nt+threadNum+Nt*threadNum+threadNum)*sizeof(cuComplex)+threadNum*sizeof(float);
				sharedMem=6000*sizeof(cuComplex);
//	cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapsize);
//	int Eu_num;
//	memset(Eu_num,0,sizeof(int));
//	Eu_num=0;
   float duration;
   clock_t start, end;
//   cudaProfilerStart();
   start=clock();
	FEpath<<<blockNum, threadNum,sharedMem>>>(d_R, d_s_hat,d_s_potential_matrix,d_s_sub_index, d_Eu, Nr, Nt, M,threadNum,d_list,d_psymbolconstellation);
	end=clock();
//	cudaProfilerStop();
	duration=double(end-start);
	printf("hey %0.4f ", duration);
	printf("\n");
	error=cudaDeviceSynchronize();
		if(error!=cudaSuccess)
		{
//		printf("cudaDeviceSynchronize returned error code %d, line %d\n", error, __LINE__);
//				 	exit(EXIT_FAILURE);
		}
		printf("%s\n",cudaGetErrorString(cudaGetLastError()));
//	printf("Eu_num is %d", Eu_num);
    error=cudaMemcpy(s_potential_matrix,d_s_potential_matrix,blockNum*Nt*sizeof(cuComplex),cudaMemcpyDeviceToHost);
    if(error!=cudaSuccess)
    {
//   	printf("s_potential_matrix returned error code %d, line %d\n", error, __LINE__);
//    			 	exit(EXIT_FAILURE);
   	}
    printf("all the potential symbol vector is:\n");
    for(count1=0;count1<blockNum;count1++)
    {
    	for(int count2=0;count2<Nt;count2++)
    	{
    		printf("%0.4f%+0.4fi ", s_potential_matrix[IDC2D(count1,count2,Nt)].x,s_potential_matrix[IDC2D(count1,count2,Nt)].y);
    	}
    	printf("\n");
    }
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

    error=cudaMemcpy(Eu,d_Eu,blockNum*sizeof(float),cudaMemcpyDeviceToHost);
    if(error!=cudaSuccess)
    {
//   	printf("Eu returned error code %d, line %d\n", error, __LINE__);
//    			 	exit(EXIT_FAILURE);
   	}


 //fine out the symbol vector index among all the block output Euclidean distance
    int Eu_mini_index=0;
    float Eu_mini_value=Eu[0];
    for(count1=0;count1<blockNum;count1++)
    {
      if(Eu[count1]<Eu_mini_value)
      {
    	  Eu_mini_value=Eu[count1];
    	  Eu_mini_index=count1;
      }
    }
    for(count1=0;count1<Nt;count1++)
    {
     s_kernel[count1]=s_potential_matrix[IDC2D(Eu_mini_index,count1,Nt)];
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
//			   	cudaFree(d_R);

}

