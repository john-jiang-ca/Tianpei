/*
 * common.h
 *
 *  Created on: Sep 12, 2014
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef COMMON_H_
#define COMMON_H_
#include <gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector_complex.h>
#include <gsl/gsl_vector_ulong.h>
#include<cuComplex.h>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#define MATRIX_SIZE 32
#define IDC2D(i,j,row)    (i)*row+j
#define Constellationsize 16
void FCSD_detection(
		cuComplex *sigRec,   //received signal vector
		cuComplex *pH,        //propagation matrix
		int Nt,             //number of transmit antennas
		int Nr,              //number of receive antennas
		int M,               //modulation scheme
		cuComplex *psymbolconstellation, //the symbol constellation
		float SNR,          //signal to noise ratio
		cuComplex *symOut,    //output symbol vector
		float *durationKernel
		);
void FCSD_decoding(
		cuComplex *d_R,  //upper triangular matrix after cholesky factorization store in device side
//		cuComplex *s_sub, //the sub brute force rho vector matrix
		cuComplex *s_hat,  //unconstrained estimation of transmitted symbol vector s
		cuComplex *s_kernel,  //quantization of estimation ,decoding results
//		cuComplex *Eu,  //Euclidean distance
		int Nt,    //the number of transmit antennas
		int Nr,    //the number of receive antennas
		int M,    //modulation scheme
		int *list,   //the permutation list
		cuComplex *psymbolconstellation //the symbol constellation
		);
void fullfact(
		int rho,  //the number of elements that use full expansion
		int M,    //constellation size
		int *s_sub_index    //the indexes of the full expansion matrix  (pow(M,rho))
		);
void comb(int m, int k, int row, int * aaa, int *subset, int *count3);
void MATRIX_INVERSE(
	cuComplex *H,  //input square matrix
	cuComplex *R,   //the inversion of the matrix H
	int Nr,  //number of receive antennas
	int Nt   //number of transmit antennas
);
void FCSD_ordering(
		cuComplex *pH,
		int *list,
		cuComplex *pH_permuted
);
void chol(cuComplex *U
);
void FCSD_CPU(
		gsl_vector_complex *preceived,
		gsl_matrix_complex *pH,
				int Nt,
				int Nr,              //number of receive antennas
				int M,               //modulation scheme
		gsl_vector_complex *psymbolconstellation, //the symbol constellation
				float SNR,
		gsl_vector_complex *symOut,
		float *durationKernel
);
void FCSD_ordering_CPU(
		cuComplex *pH,
		int *list,
		cuComplex *pH_permuted
		);


__global__ void queue_max(
		cuComplex *P,
		int *j,
		int index
		);
__global__ void queue_min(
		cuComplex *P,
		int *j,
		int index
		);
__global__ void MATRIX_ROWCOLUMNT_kernel(cuComplex *Hr,
		cuComplex *Hc, int row, int column);
__global__ void MATRIX_COLUMNROWT_kernel(cuComplex *Hc,
		cuComplex *Hr, int row, int column);
#endif /* COMMON_H_ */
