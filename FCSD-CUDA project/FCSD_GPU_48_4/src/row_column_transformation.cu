/*
 * row_column_transformation.cu
 *
 *  Created on: Jul 7, 2014
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include"common.h"
#include<cuComplex.h>
#include<cuda.h>
#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include"cu_complex_operation.cuh"



__global__ void MATRIX_ROWCOLUMNT_kernel(cuComplex *Hr,
		cuComplex *Hc, int row, int column)
{
	int i=blockIdx.x;
	int j=threadIdx.x;


     Hc[IDC2D(j,i,row)]=Hr[IDC2D(j,i,column)];
     __syncthreads();

}




__global__ void MATRIX_COLUMNROWT_kernel(cuComplex *Hc,
		cuComplex *Hr, int row, int column)
{
	int i=blockIdx.x;
	int j=threadIdx.x;


     Hr[IDC2D(i,j,column)]=Hc[IDC2D(i,j,row)];
     __syncthreads();


}





void MATRIX_COLUMNROWT(cuComplex *Hc, cuComplex *Hr, int row, int column){

	cuComplex *d_pCC, *d_pCCR;

	cudaMalloc((void**)&d_pCCR,row*column*sizeof(cuComplex));
	cudaMalloc((void**)&d_pCC,row*column*sizeof(cuComplex));
	cudaMemcpy(d_pCC,Hc,row*column*sizeof(cuComplex),cudaMemcpyHostToDevice);
	dim3 blockDim(column,1);
	MATRIX_COLUMNROWT_kernel<<<row,column,0>>>(d_pCC, d_pCCR, row, column);
	cudaDeviceSynchronize();
	cudaMemcpy(Hr,d_pCCR,row*column*sizeof(cuComplex),cudaMemcpyDeviceToHost);
	cudaFree(d_pCC);
	cudaFree(d_pCCR);

}
void MATRIX_ROWCOLUMNT(cuComplex *Hr, cuComplex *Hc, int row,int column){

	cuComplex *d_pCC, *d_pCCR;

	cudaMalloc((void**)&d_pCCR,row*column*sizeof(cuComplex));
	cudaMalloc((void**)&d_pCC,row*column*sizeof(cuComplex));
	cudaMemcpy(d_pCC,Hr,row*column*sizeof(cuComplex),cudaMemcpyHostToDevice);
	dim3 blockDim(column,1);
	MATRIX_ROWCOLUMNT_kernel<<<row,column,0>>>(d_pCC, d_pCCR, row,column);
	cudaDeviceSynchronize();
	cudaMemcpy(Hc,d_pCCR,row*column*sizeof(cuComplex),cudaMemcpyDeviceToHost);
	cudaFree(d_pCC);
	cudaFree(d_pCCR);

}
