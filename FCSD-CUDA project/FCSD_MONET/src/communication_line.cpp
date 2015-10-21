/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
/**
 *This is the communication line performing the fixed complexity sphere decoding algorithm
 * * This function implement a complete communication line including: data generation, modulation,
 * channel generation, detection, decoding, print the result, speed and bit error rate
 * Tianpei chen
 * Email:tianpei.chen@mail.mcgill.ca
 * 2014.08.19
 * pdata: the transmitted gray data in decimal
 * pgraycode: find the indexes of symbol according to graycode
 * pgrayindex: find the graycode according to the indexes of the symbol
 * pgraydata: find the indexes according to the pdata
 * ptransmitted: find the complex symbol corresponding to he pgraydata
 */
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_ulong.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
/*#include <gsl/gsl_matrix.h>*/
#include <gsl/gsl_matrix_complex_float.h>
#include <gsl/gsl_randist.h>
/*#include <gsl/gsl_cblas.h>*/
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
//#include<cuda_runtime.h>
#include<assert.h>
#include<time.h>
#include"channel_generator.h"
#include"grayencoder.h"
#include"symbolconstellation.h"
#include"modulator.h"
#include"noise_generator.h"
//#include"chol.cuh"
#include"data_generator.h"
//#include"demodulator.h"
#include"binaryerrors.h"
#include"complex.h"
//#include"chol.cuh"
//#include"cu_complex_operation.cuh"
//#include "helper_functions.h"
//#include"chol_common.h"
#include<time.h>
//#include<cublas_v2.h>
//#include<cuComplex.h>
//#include<cuda_runtime.h>
//#include<cuda.h>
//#include"cu_complex_operation.cuh"
//#include<cudaProfiler.h>
//#include<cuda_profiler_api.h>
#include "common.h"
#define BER "bit_error_rate.txt"   //the file to store the bit error rate
#define SER "symbol_error_rate.txt"  //the file to store the symbol error rate
#define timeused "time.txt"              //the file to store operation time of detection algorithm
#define GFLOPS "GFLOPS.txt"          //the file to store the Giga flops operated per second
#define miniteration 1e1               //the minimum number of channel realizations
#define minSymbolError 50          //the minimum number of symbol error
#define epsilon 1e-5               //the accuracy
void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, unsigned long Q);
void grayencoder(gsl_vector_ulong *pgraydata, gsl_vector_ulong *pgrayindexes,
		unsigned long Q);
void symbolconstellation(gsl_vector_complex *psymbolconstellation, int Nt);
void modulator(gsl_vector_ulong *pdata, gsl_vector_complex *ptransmitted,
		gsl_vector_complex *psymbolconstellation);
void channel_generator(gsl_matrix_complex *pH, gsl_rng *pr);
//void error_channel_generator (gsl_matrix_complex *pH, gsl_rng *pr, double sigman);
void noise_generator(gsl_vector_complex *pnoise, gsl_rng *pr, double sigman);
//void demodulator(cuComplex *symOut, gsl_vector_complex *psymbolconstellation,gsl_vector_ulong *pgrayindexes,gsl_vector_ulong *poutput);
//void matrix_invert_complex (gsl_matrix_complex *start, gsl_matrix_complex *result);
//void matrix_sqrt (gsl_matrix_complex *pstart, gsl_matrix_complex *presult);
//void corr_matrices (gsl_matrix_complex *pRt, gsl_matrix_complex *pRr, double Rt, double Rr);
unsigned long binaryerrors(unsigned long input, unsigned long output,
		unsigned long Q);
//void GreedyMMSE (gsl_matrix_complex *pH, gsl_vector_complex *preceived, gsl_vector_int *porder, gsl_matrix_complex *pR, gsl_vector_complex *pytilde, gsl_matrix_complex *pHord, double alpha, int N, double *pPreFlops);
//void FSDtrianginv (gsl_matrix_complex *pH, gsl_matrix_complex *pR, gsl_vector_int *porder, int N, double *pPreFlops);
//void Givensupdate (gsl_matrix_complex *pR, gsl_vector_complex *pytilde, double *pPreFlops);
//void rings (gsl_vector_complex *psymbolconstellation, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int *numberofrings);
//void constrcandidates (gsl_comple256x sbabai, double bound, double d, int M, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int numberofrings, gsl_vector_complex *psymbolconstellation,  int maxcandidates, double alpha, gsl_complex rll, int *pfoundcandidates, gsl_vector_complex *pcandidates, gsl_vector_ulong *poutputcandidates, gsl_vector *pdistances, int *ptotalcandidates);
//void Kbest (gsl_vector_complex *pytilde, gsl_matrix_complex *pR, gsl_vector_int *pchildren, int N, gsl_vector_int *psurvivors, double alpha, double initradiussq, gsl_vector_complex *psymbolconstellation, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int numberofrings, gsl_vector_ulong *poutput);

int main(void) {
	int iteration;                       //iteration time
	int bitError, symError;               //bit error symbol error
	float BitErrorRate, SymbolErrorRate, sigmas, sigman; //bit error rate, symbol error rate, standard deviation of signal and noise
	float *SNR;
	int M = Constellationsize;                               //modulation scheme
	int count, count1, count2, count3, count4;   //used for loops
	const int Nt = MATRIX_SIZE; // Nt is the number of transmit antennas. Must be a positive integer.
	const int Nr = MATRIX_SIZE; // Nr is the number of receive antennas. Must be a positive integer.
	unsigned long temp, errors, symbolError;
	symbolError = 0;
    gsl_vector_complex *symOut;
	symOut = gsl_vector_complex_calloc(Nt);
	gsl_complex one;
	GSL_SET_COMPLEX(&one, 1.0, 0.0256);
	gsl_complex zero;
	GSL_SET_COMPLEX(&zero, 0.0, 0.0);
	SNR = (float*) malloc(8 * sizeof(float));
	for (count1 = 0; count1 < 8; count1++) {
		SNR[count1] = pow(10, ((count1 + 1) * 2 / 10)); //test the SNR from 2 to 16 the step is 2
	}

	//	switch(M)
	//	{
	//	case 2;     //BPSKs_potential_matrix[IDC2D(blockIdx.x,(list[threadIdx.x]-1),Nt)]
	//	case 4;     //QPSK
	//	case 16;   //16QAM
	//	case 64;   //64QAMount1=0;count1<MAT
	//	}
	FILE *cfile1, *cfile2, *cfile3, *cfile4;
	cfile1 = fopen(BER, "a");
	cfile2 = fopen(SER, "a");
	cfile3 = fopen(timeused, "a");
	cfile4 = fopen(GFLOPS, "a");
	//gray encoding and constellation definition
	gsl_vector_ulong *pgraycode, *pgrayindexes;
	gsl_vector_complex *psymbolconstellation;
	gsl_rng *pr;
	pgraycode = gsl_vector_ulong_calloc(M);
	pgrayindexes = gsl_vector_ulong_calloc(M);
	psymbolconstellation = gsl_vector_complex_calloc(M);
	grayencoder(pgraycode, pgrayindexes, M);
	symbolconstellation(psymbolconstellation, Nt);
	unsigned long seed = 0L;
	const gsl_rng_type *pT;
	pT = gsl_rng_default;
	pr = gsl_rng_alloc(pT);
	seed = time(NULL);
	gsl_rng_set(pr, seed);
	//	for(count1=0;count1<MATRIX_SIZE;count1++)
	clock_t time1, time2;
	double duration_total;
	float *durationKernel_CPU, duration_CPU;
	float *durationKernel_GPU, duration_GPU;
	float durationKernel_CPU_t = 0;
	float durationKernel_GPU_t = 0;
	durationKernel_CPU = (float*) malloc(sizeof(float));
	durationKernel_GPU = (float*) malloc(sizeof(float));
	time1 = clock();
	for (count = 0; count < 100; count++) {
		iteration = 0;
		bitError = 0;
		symError = 0;
		BitErrorRate = 0;
		SymbolErrorRate = 0;
		//	do
		//	{

		gsl_vector_ulong *pdata, *pgraydata, *poutput;
		gsl_vector_complex *ptransmitted, *pnoise, *preceived, *symOut;
		gsl_matrix_complex *pH;
		//		cuComplex *ptransmitted_cu,*psymbolconstellation_cu,*sigRec,*pH_cu,*symOut_cu;
		pdata = gsl_vector_ulong_calloc(Nt);
		pgraydata = gsl_vector_ulong_calloc(Nt);
		poutput = gsl_vector_ulong_calloc(Nt);
		ptransmitted = gsl_vector_complex_calloc(Nt);
		pnoise = gsl_vector_complex_calloc(Nr);
		preceived = gsl_vector_complex_calloc(Nr);
		pH = gsl_matrix_complex_calloc(Nr, Nt);
		symOut = gsl_vector_complex_calloc(Nt);
		//		ptransmitted_cu=(cuComplex*)malloc(Nt*sizeof(cuComplex));
		//		psymbolconstellation_cu=(cuComplex*)malloc(M*sizeof(cuComplex));
		//		sigRec=(cuComplex*)malloc(Nr*sizeof(cuComplex));
		//		pH_cu=(cuComplex*)malloc(Nr*Nt*sizeof(cuComplex));
		//		symOut_cu=(cuComplex*)malloc(Nt*sizeof(cuComplex));
		//data generation and channel generation
		data_generator(pdata, pr, M); //pdata is actually gray code
		for (count2 = 0; count2 < Nt; count2++) {
			temp = gsl_vector_ulong_get(pdata, count2);
			gsl_vector_ulong_set(pgraydata, count2,
					gsl_vector_ulong_get(pgraycode, temp)); //pgraydata is the indexes corresponding to gray code
			//			gsl_vector_ulong_set(pgraydata,count2,gsl_vector_ulong_get(pdata,count2));   //find the gray
		}
		modulator(pgraydata, ptransmitted, psymbolconstellation); //ptransmitted stores the complex symbol

		channel_generator(pH, pr);
		sigmas = sqrt(1.0 / Nt);
		sigman = sqrt(pow(sigmas, 2) / (SNR[7])); /* corresponding noise standard deviation per dimension */
		noise_generator(pnoise, pr, sigman);
		for (count1 = 0; count2 < Nt; count2++) {
			gsl_vector_complex_set(preceived, count1,
					gsl_vector_complex_get(pnoise, count1));
		}
		gsl_blas_zgemv(CblasNoTrans, one, pH, ptransmitted, one, preceived);
		//		for(count2=0;count2<M;count2++)
		//		{
		//			psymbolconstellation_cu[count2].x=gsl_vector_complex_get(psymbolconstellation,count2).dat[0];
		//			psymbolconstellation_cu[count2].y=gsl_vector_complex_get(psymbolconstellation,count2).dat[1];
		//		}
		//		for( count1=0;count1<M;count1++)
		//				{
		//					printf("the symbol constellation is: %0.4f%+0.4fi ",psymbolconstellation_cu[count1].x,psymbolconstellation_cu[count1].y );
		//				}
		//		printf("\n");
		//		for(count2=0;count2<Nt;count2++)
		//		{
		//			ptransmitted_cu[count2].x=gsl_vector_complex_get(ptransmitted,count2).dat[0];
		//			ptransmitted_cu[count2].y=gsl_vector_complex_get(ptransmitted,count2).dat[1];
		//		}
		//		for(count2=0;count2<Nr;count2++)
		//		{
		//			sigRec[count2].x=gsl_vector_complex_get(preceived,count2).dat[0];
		//			sigRec[count2].y=gsl_vector_complex_get(preceived,count2).dat[1];
		//			for(count3=0;count3<Nt;count3++)
		//			{
		//				pH_cu[IDC2D(count2,count3,Nt)].x=gsl_matrix_complex_get(pH,count2,count3).dat[0];
		//				pH_cu[IDC2D(count2,count3,Nt)].y=gsl_matrix_complex_get(pH,count2,count3).dat[1];
		//			}
		//		}
		double start, end, duration1, duration2;
		start = clock();
		FCSD_CPU(preceived, pH, Nt, Nr, M, psymbolconstellation, SNR[7], symOut,
				durationKernel_CPU);
		end = clock();
		duration1 = double(end - start);
		durationKernel_CPU_t += *durationKernel_CPU;
		duration_CPU += duration1;

		//		  cudaProfilerStart();
		start = clock();
		//        FCSD_detection(sigRec, pH_cu,  Nt, Nr,  M, psymbolconstellation_cu, SNR[7],  symOut_cu, durationKernel_GPU );   //detection algorithm
		end = clock();

		duration2 = double(end - start);
		durationKernel_GPU_t += *durationKernel_GPU;
		//         cudaProfilerStop();
		duration_GPU += duration2;
		printf("the operation time in GPU is:\n");
		printf("%0.4f ", duration2);
		printf("the operation time in CPU is:\n");
		printf("%0.4f ", duration1);
		printf("the original transmitted symbol vector is:\n");
		//        for(count1=0;count1<MATRIX_SIZE;count1++)
		//        {
		//        	printf("%0.4f%+0.4fi ", ptransmitted_cu[count1].x, ptransmitted_cu[count1].y);
		//        }
		printf("the decoded transmitted symbol vector by CPU is:\n");
		for (count1 = 0; count1 < MATRIX_SIZE; count1++) {
			printf("%0.4f%+0.4fi ",
					gsl_vector_complex_get(symOut, count1).dat[0],
					gsl_vector_complex_get(symOut, count1).dat[1]);
		}
		printf("the decoded transmitted symbol vector by GPU is:\n");
		//        for(count1=0;count1<MATRIX_SIZE;count1++)
		//        {
		//        	printf("%0.4f%+0.4fi ", symOut_cu[count1].x, symOut_cu[count1].y);
		//        }
		//        demodulator(symOut_cu, ptransmitted,pgrayindexes,poutput);  //poutput is the graycode of the symbols
		iteration = iteration + 1;
		//check symbol error
		//         for(count3=0;count3<Nt;count3++)
		//         {
		//        	 if(abs(symOut_cu[count3].x-ptransmitted_cu[count3].x)>epsilon||abs(symOut_cu[count3].y-ptransmitted_cu[count3].y)>epsilon)
		//        	 {
		//        		 symbolError+=symbolError;
		//        	 }
		//         }
		//check bit error
		for (count4 = 0; count4 < Nt; count4++) {
			//        	 errors = binaryerrors(gsl_vector_ulong_get (pdata, count4), gsl_vector_ulong_get (poutput,count4), M);
			errors = 1;
			bitError += errors;
		}
		printf("the iteration time is %d\n", count1);
		gsl_vector_ulong_free(pdata);
		pdata = NULL;
		gsl_vector_ulong_free(pgraydata);
		pgraydata = NULL;
		gsl_vector_ulong_free(poutput);
		poutput = NULL;
		gsl_vector_complex_free(ptransmitted);
		ptransmitted = NULL;
		gsl_vector_complex_free(pnoise);
		pnoise = NULL;
		gsl_vector_complex_free(preceived);
		preceived = NULL;
		gsl_matrix_complex_free(pH);
		pH = NULL;
		gsl_vector_complex_free(symOut);
		symOut = NULL;
		//         	free(ptransmitted_cu);
		//         	ptransmitted_cu=NULL;
		//         	free(psymbolconstellation_cu);
		//         	psymbolconstellation_cu=NULL;
		//         	free(sigRec);
		//         	sigRec=NULL;
		//         	free(pH_cu);
		//         	pH_cu=NULL;
		//         	free(symOut_cu);
		//         	symOut_cu=NULL;
		//	}
		//     while ((symError < minSymbolError)||(iteration < miniteration));

	}
	time2 = clock();
	//	printf("the clock per second is %f:\n", CLOCKS_PER_SEC);
	duration_total = double(time2 - time1) / CLOCKS_PER_SEC;
	fprintf(cfile3, "this is for monet 8X8\n");
	fprintf(cfile3, "the total time is= %f \n", duration_total);
	fprintf(cfile3, "the total duration of CPU kernel is %f: \n",
			double(durationKernel_CPU_t / double(CLOCKS_PER_SEC)));
	fprintf(cfile3, "the total duration of GPU kernel is %f: \n",
			double(durationKernel_GPU_t / double(CLOCKS_PER_SEC)));
	fprintf(cfile3, "the total duration of CPU is %f: \n",
			double(duration_CPU / double(CLOCKS_PER_SEC)));
	fprintf(cfile3, "the total duration of GPU is %f: \n",
			double(duration_GPU / double(CLOCKS_PER_SEC)));
	gsl_rng_free(pr);
	pr = NULL;
	gsl_vector_ulong_free(pgraycode);
	pgraycode = NULL;
	gsl_vector_ulong_free(pgrayindexes);
	pgrayindexes = NULL;
	gsl_vector_complex_free(psymbolconstellation);
	psymbolconstellation = NULL;
	free(SNR);
	free(durationKernel_CPU);
	free(durationKernel_GPU);
	SNR = NULL;
	printf("the whole program finishes successfully\n");
	fclose(cfile1);
	fclose(cfile2);
	fclose(cfile3);
	fclose(cfile4);
	return (0);

}
