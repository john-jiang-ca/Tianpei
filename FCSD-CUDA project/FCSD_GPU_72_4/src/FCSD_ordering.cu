/*
 * FCSD_ordering.c
 *this function implement the channel ordering of fixed complexity sphere decoding algorithm
 *INPUT:
 *pH: original propagation channel matrix
 *list: permutation list
 *pH_permuted: permuted propagation channel matrix
 *  Created on: Jun 26, 2014
 *      Author: Tianpei Chen
 *      Email: tianpei.chen@mail.mcgill.ca
 */

//#ifndef FCSD_ORDERING_H_
//#define FCSD_ORDERING_H_
#include"common.h"
#include<cublas_v2.h>
#include<cuComplex.h>
#include<cuda_runtime.h>
#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<assert.h>
#include<string.h>
#include<cuda.h>
void FCSD_ordering(
		cuComplex *pH,
		int *list,
		cuComplex *d_pH_permuted
)
{

	int count, count1,count2,count3;
	int mm,*j,*list_temp;
	j=(int*)malloc(sizeof(int));
	list_temp=(int*)malloc(MATRIX_SIZE*sizeof(int));
//	int *d_j;
//	cudaMalloc((void**) &d_j, sizeof(int));
	cuComplex *pH_loop;
	pH_loop=(cuComplex*)malloc(MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
	int N=ceil(sqrt(MATRIX_SIZE)-1);//
	cublasHandle_t handle;
	cublasStatus_t ret;
	cudaError_t error;
	ret=cublasCreate(&handle);

    cuComplex *d_pH, *d_Pprod, *d_Pinv;
    cuComplex alpha,beta;
    alpha.x=1;
    alpha.y=0;
    beta.x=0;
    beta.y=0;
	  clock_t start,end;
	  start=clock();
		error=cudaMalloc((void**) &d_Pprod, (MATRIX_SIZE)*(MATRIX_SIZE)*sizeof(cuComplex));//store the wishart matrix
			if(error!=cudaSuccess)
			{
				printf("cudaMalloc d_Pprod returned error code %d, line %d\n", error, __LINE__);
//				exit(EXIT_FAILURE);
			}
			error=cudaMalloc((void**) &d_Pinv, (MATRIX_SIZE)*(MATRIX_SIZE)*sizeof(cuComplex));//store the inverse of wishart matrix
			if(error!=cudaSuccess)
			{
				printf("cudaMalloc d_Pinv returned error code %d, line %d\n", error, __LINE__);
				exit(EXIT_FAILURE);
			}
			cudaMalloc((void**) &d_pH,(MATRIX_SIZE)*MATRIX_SIZE*sizeof(cuComplex));



	for(count2=0;count2<MATRIX_SIZE-1;count2++)
	{
			if(count2==0)
			{
			cudaMemcpy(d_pH,pH,(MATRIX_SIZE-count2)*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyHostToDevice);
			}
			else
			{
			cudaMemcpy(d_pH,pH_loop,(MATRIX_SIZE-count2)*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyHostToDevice);
			}


		ret=cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_C, MATRIX_SIZE-count2, MATRIX_SIZE-count2, MATRIX_SIZE, &alpha, d_pH, MATRIX_SIZE-count2, d_pH, MATRIX_SIZE-count2, &beta, d_Pprod, MATRIX_SIZE-count2);
	    if (ret != CUBLAS_STATUS_SUCCESS)
	    {
	        printf("cublasSgemm returned error code %d, line(%d)\n", ret, __LINE__);
	        exit(EXIT_FAILURE);
	    }



	    int column=MATRIX_SIZE-count2;
	    int row=MATRIX_SIZE;
	    double duration;

	    start=clock();
	    MATRIX_INVERSE(d_Pprod,d_Pinv,row,column);
	    end=clock();
	    duration=double((end-start));
	    cudaDeviceSynchronize();
	    printf("the duration of matrix inverse GPU is %0.4f:\n", duration);
		if(count2<N)
		{
			cublasIcamax(handle, (MATRIX_SIZE-count2),d_Pinv,(MATRIX_SIZE-count2),j);
		}
		else
		{
			cublasIcamin(handle, (MATRIX_SIZE-count2),d_Pinv,(MATRIX_SIZE-count2),j);
		}

		  list_temp[count2]=list[count2+*j-1];   //choose the *j th element in the new list to be detected first
					    list[count2+*j-1]=0;
					for (count3=0;count3<MATRIX_SIZE-count2-1;count3++)
					{//				double durationMalloc;
						//				start=clock();
						//				cudaMalloc((void**) &d_pH,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
						//				cudaMemcpy(d_pH,pH,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyHostToDevice);
						////				cudaMemcpy(pH,d_pH,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyDeviceToHost);
						//			    end=clock();
						//			    durationMalloc=double(end-start);
						//			    printf("the memcpy time is %0.4f:\n", durationMalloc);



						//		error=cudaMalloc((void**) &d_pH, (MATRIX_SIZE-count2)*MATRIX_SIZE*sizeof(cuComplex));// store the Hi after each refresh
						//		if(error!=cudaSuccess)
						//		{
						//			printf("cudaMalloc d_pH returned error code %d, line %d\n", error, __LINE__);
						//			exit(EXIT_FAILURE);
						//		}
						for(int count4=count2;count4<MATRIX_SIZE;count4++)
						{
							if(list[count4]!=0)
							{
						list_temp[count2+count3+1]=list[count4];   //choose the next MATRIX_SIZE-(count2+1) to form the next sequence, the sequence order is the same
						list[count4]=0;
							break;
							}
						}
					}
					for(count3=0;count3<MATRIX_SIZE;count3++)
					{
					list[count3]=list_temp[count3];
					}

	cuComplex *pH_temp1=(cuComplex*)malloc((MATRIX_SIZE)*(MATRIX_SIZE-count2-1)*sizeof(cuComplex));

						for(count3=0;count3<MATRIX_SIZE;count3++)
						 {
							 for(count1=0;count1<MATRIX_SIZE-count2;count1++)
							 {
								 if(count1<(*j-1))
								 {
			//		              pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]=pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
									 pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]=pH[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
								 }
								 else if(count1>(*j-1))
								 {
			//						 pH_temp1[IDC2D(count3,count1-1,(MATRIX_SIZE-count2-1))]=pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
									 pH_temp1[IDC2D(count3,(count1-1),(MATRIX_SIZE-count2-1))]=pH[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
								 }
							 }
						 }
pH_loop=(cuComplex*)realloc(pH_loop,(MATRIX_SIZE)*(MATRIX_SIZE-count2-1)*sizeof(cuComplex));
						for(count3=0;count3<MATRIX_SIZE;count3++)
						 {
							 for(count1=0;count1<MATRIX_SIZE-count2-1;count1++)
							 {
			//		         pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]= pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))];
								pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]=pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))];
							 }
						 }


free(pH_temp1);
}
end=clock();
double duration=double(end-start);
	cuComplex *pH_permuted;
	pH_permuted=(cuComplex*)malloc(MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
	for(count1=0;count1<MATRIX_SIZE;count1++)
	{
	for(count2=0;count2<MATRIX_SIZE;count2++)
	{
	pH_permuted[IDC2D(count1,(MATRIX_SIZE-count2-1),MATRIX_SIZE)]=pH[IDC2D(count1,list[count2]-1,MATRIX_SIZE)];
	}
	}
  cudaMemcpy(d_pH_permuted,pH_permuted,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyHostToDevice);
    ret=cublasDestroy(handle);
    if (ret != CUBLAS_STATUS_SUCCESS)
       {
           printf("cublasDestory returned error code %d, line(%d)\n", ret, __LINE__);
           exit(EXIT_FAILURE);
       }

cudaFree(d_pH);
cudaFree(d_Pprod);
cudaFree(d_Pinv);
free(list_temp);
free(j);
free(pH_loop);
free(pH_permuted);

}



