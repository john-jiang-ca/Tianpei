//============================================================================
// Name        : SVR_MIMO_communication_line.cpp
// Author      : Tianpei Chen
// Version     :
// Copyright   : Your copyright notice
// Description : This is the main routine for large MIMO communication line,
// including random data generator, gray encoder, modulator, channel generator, noise generator
// detector, symbol error and binary error accumulator, SER and BER calculator.
//============================================================================

#include"channel_generator.h"
#include"grayencoder.h"
#include"symbolconstellation.h"
#include"modulator.h"
#include"noise_generator.h"
#include"data_generator.h"
#include"demodulator.h"
#include"binaryerrors.h"
#include"SVR_DETECTOR.h"
#include "public.h"
#define BER_SVR "BER.txt"   //the file to store the bit error rate
#define SER_SVR "SER.txt"  //the file to store the symbol error rate
#define TIME_SVR "time.txt"              //the file to store operation time of detection algorithm
#define GFLOPS_SVR "GFLOPS.txt"          //the file to store the Giga flops operated per second
#define miniteration 100               //the minimum number of channel realizations
#define minSymbolError 5          //the minimum number of symbol error
#define SNRnum  10                    //the point number of signal to noise ratio per bit



void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, unsigned long Q);
void grayencoder (gsl_vector_ulong *pgraydata, gsl_vector_ulong *pgrayindexes, unsigned long Q);
void symbolconstellation(gsl_vector_complex *psymbolconstellation, int Nt_conste);
void modulator (gsl_vector_ulong *pdata, gsl_vector_complex *ptransmitted,
		gsl_vector_complex *psymbolconstellation);
void channel_generator (gsl_matrix_complex *pH, gsl_rng *pr);
void noise_generator (gsl_vector_complex *pnoise, gsl_rng *pr, double sigman);
void SVR_DETECTOR( gsl_vector_complex *preceived,  //received symbol vector
	 gsl_matrix_complex *pH, //channel propagation matrix
		float SNRb, // bit signal to noise ratio
		gsl_vector_complex *psymout, //detected transmitted symbol vector
		int start_M,  //initialization method
		int selection_M //work set selection strategy
		);
void demodulator_CPU(gsl_vector_complex *symOut, gsl_vector_complex *psymbolconstellation,
		gsl_vector_ulong *pgrayindexes, gsl_vector_ulong *poutput);
unsigned long binaryerrors (unsigned long input, unsigned long output, unsigned long Q);
//void error_channel_generator (gsl_matrix_complex *pH, gsl_rng *pr, double sigman);

//void demodulator_GPU(cuComplex *symOut, gsl_vector_complex *psymbolconstellation,
//		gsl_vector_ulong *pgrayindexes, gsl_vector_ulong *poutput);
//void matrix_invert_complex (gsl_matrix_complex *start, gsl_matrix_complex *result);
//void matrix_sqrt (gsl_matrix_complex *pstart, gsl_matrix_complex *presult);
//void corr_matrices (gsl_matrix_complex *pRt, gsl_matrix_complex *pRr, double Rt, double Rr);

//void GreedyMMSE (gsl_matrix_complex *pH, gsl_vector_complex *preceived, gsl_vector_int *porder, gsl_matrix_complex *pR, gsl_vector_complex *pytilde, gsl_matrix_complex *pHord, double alpha, int N, double *pPreFlops);
//void FSDtrianginv (gsl_matrix_complex *pH, gsl_matrix_complex *pR, gsl_vector_int *porder, int N, double *pPreFlops);
//void Givensupdate (gsl_matrix_complex *pR, gsl_vector_complex *pytilde, double *pPreFlops);
//void rings (gsl_vector_complex *psymbolconstellation, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int *numberofrings);
//void constrcandidates (gsl_comple256x sbabai, double bound, double d, int M, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int numberofrings, gsl_vector_complex *psymbolconstellation,  int maxcandidates, double alpha, gsl_complex rll, int *pfoundcandidates, gsl_vector_complex *pcandidates, gsl_vector_ulong *poutputcandidates, gsl_vector *pdistances, int *ptotalcandidates);
//void Kbest (gsl_vector_complex *pytilde, gsl_matrix_complex *pR, gsl_vector_int *pchildren, int N, gsl_vector_int *psurvivors, double alpha, double initradiussq, gsl_vector_complex *psymbolconstellation, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int numberofrings, gsl_vector_ulong *poutput);

int main(void) {
	long int iteration;                       //iteration time
	double BitErrorRate, SymbolErrorRate, sigman; //bit error rate, symbol error rate, standard deviation of signal and noise
	float *SNR;      //this is the SNR per bit
	int count, count1, count2, count3, count4;   //used for loops
//	const int Nt = MATRIX_SIZE; // Nt is the number of transmit antennas. Must be a positive integer.
//	const int Nr = MATRIX_SIZE; // Nr is the number of receive antennas. Must be a positive integer.
	unsigned long temp, symbolError, BitError;
	gsl_complex one;
	GSL_SET_COMPLEX(&one, 1.0, 0.0);
	gsl_complex zero;
	GSL_SET_COMPLEX(&zero, 0.0, 0.0);
	SNR = (float*) malloc(11 * sizeof(float));
	for (count1 = 0; count1 <=SNRnum; count1++) {
		SNR[count1] = pow(10, (float)(((count1) * 2) / (float)(10))); //test the SNR from 2 to 20 the step is 2
	}
	float begin=0;      //the SNR per bit in dB beginer begin=SNRnum for single SNR point begin=0 for multiple SNR points
	FILE *cfile1, *cfile2, *cfile3, *cfile4;
	cfile1 = fopen(BER_SVR, "a");
	cfile2 = fopen(SER_SVR, "a");
	cfile3 = fopen(TIME_SVR, "a");
	cfile4 = fopen(GFLOPS_SVR, "a");
	//gray encoding and constellation definition
	gsl_vector_ulong *pgraycode, *pgrayindexes;
	gsl_vector_complex *psymbolconstellation;
	gsl_rng *pr;
	pgraycode = gsl_vector_ulong_calloc(M);
	pgrayindexes = gsl_vector_ulong_calloc(M);
	psymbolconstellation = gsl_vector_complex_calloc(M);
	grayencoder(pgraycode, pgrayindexes, M);
	symbolconstellation(psymbolconstellation, Nt);  //generate the symbol constellation based on indexes
	unsigned long seed = 0L;
	const gsl_rng_type *pT;
	pT = gsl_rng_default;
	pr = gsl_rng_alloc(pT);
	seed = time(NULL);
	gsl_rng_set(pr, seed);
	clock_t start, end;
	clock_t time1, time2;
//Time cost counter in detector
	float Time_CSVR;
	float Time_SVM;
	float Time_FCSD;
	float Time_CSVR_aver[SNRnum];
	float Time_SVM_aver[SNRnum];
	float Time_FCSD_aver[SNRnum];
	float duration_total[SNRnum];
	printf("the program for %dX%d %d QAM begin!!\n", Nr, Nt,M);
	printf("the SNR per bit is:\n");
//	fprintf(cfile1,"this is the SNR %d %f\n", 2*SNRnum, SNR[SNRnum] );
	for(count1=begin;count1<=SNRnum;count1++)
	{
	printf("%d ", 2*count1);
	}
	printf("\n");
//preface
	fprintf(cfile1,"\n");
	fprintf(cfile1, "bit error rate of CSVR algorithm \n");
	fprintf(cfile1, "%d X %d MIMO\n", Nr, Nt);
	fprintf(cfile1,"%d QAM modulation\n", M);
	fprintf(cfile1, "the SNR per bit is:\n");
//	fprintf(cfile1,"this is the SNR %d %f\n", 2*SNRnum, SNR[SNRnum] );
	for(count1=begin;count1<=SNRnum;count1++)
	{
	fprintf(cfile1,"%d ", 2*count1);
	}
	fprintf(cfile1,"\n");
	for(count1=begin;count1<=SNRnum;count1++)
	{
	fprintf(cfile1,"%f ", SNR[count1]);
	}
	fprintf(cfile1,"\n");
//preface
	fprintf(cfile2,"\n");
	fprintf(cfile2, "symbol error rate of CSVR algorithm \n");
	fprintf(cfile2, "%d X %d MIMO\n", Nr, Nt);
	fprintf(cfile2,"%d QAM modulation\n", M);
	fprintf(cfile2, "the SNR per bit is:\n");
//	fprintf(cfile1,"this is the SNR %d %f\n", 2*SNRnum, SNR[SNRnum] );
	for(count1=begin;count1<=SNRnum;count1++)
	{
	fprintf(cfile2,"%d ", 2*count1);
	}
	fprintf(cfile2,"\n");
	for(count1=begin;count1<=SNRnum;count1++)
	{
	fprintf(cfile2,"%f ", SNR[count1]);
	}
	fprintf(cfile2,"\n");
//preface
		fprintf(cfile3,"\n");
		fprintf(cfile3, "Average time used for detectors \n");
		fprintf(cfile3, "%d X %d MIMO\n", Nr, Nt);
		fprintf(cfile3,"%d QAM modulation\n", M);
		fprintf(cfile3, "the SNR per bit is:\n");
	//	fprintf(cfile1,"this is the SNR %d %f\n", 2*SNRnum, SNR[SNRnum] );
		for(count1=begin;count1<=SNRnum;count1++)
		{
		fprintf(cfile3,"%d ", 2*count1);
		}
		fprintf(cfile3,"\n");
		for(count1=begin;count1<=SNRnum;count1++)
		{
		fprintf(cfile3,"%f ", SNR[count1]);
		}
		fprintf(cfile3,"\n");
//perface
		fprintf(cfile4,"\n");
		fprintf(cfile4, "Average computational cost and iteration time of detectors \n");
		fprintf(cfile4, "%d X %d MIMO\n", Nr, Nt);
		fprintf(cfile4,"%d QAM modulation\n", M);
		fprintf(cfile4, "the SNR per bit is:\n");
	//	fprintf(cfile1,"this is the SNR %d %f\n", 2*SNRnum, SNR[SNRnum] );
		for(count1=begin;count1<=SNRnum;count1++)
		{
		fprintf(cfile4,"%d ", 2*count1);
		}
		fprintf(cfile4,"\n");
		for(count1=begin;count1<=SNRnum;count1++)
		{
		fprintf(cfile4,"%f ", SNR[count1]);
		}
		fprintf(cfile4,"\n");

	for (count = begin; count <=SNRnum; count++) {
		iteration = 0;
		BitError = 0;
		symbolError = 0;
		BitErrorRate = 0;
		SymbolErrorRate = 0;
		Time_CSVR=0;
		Time_FCSD=0;
		Time_SVM=0;
		time1=clock();
		do {

			gsl_vector_ulong *pdata, *pgraydata, *poutput;
			gsl_vector_complex *ptransmitted, *pnoise, *preceived, *symOut;
			gsl_matrix_complex *pH;
			pdata = gsl_vector_ulong_calloc(Nt);
			pgraydata = gsl_vector_ulong_calloc(Nt);
			poutput = gsl_vector_ulong_calloc(Nt);
			ptransmitted = gsl_vector_complex_calloc(Nt);
			pnoise = gsl_vector_complex_calloc(Nr);
			preceived = gsl_vector_complex_calloc(Nr);
			pH = gsl_matrix_complex_calloc(Nr, Nt);
			symOut = gsl_vector_complex_calloc(Nt);
			//data generation and channel generation
			data_generator(pdata, pr, M); //pdata is the indexes
			for (count2 = 0; count2 < Nt; count2++) {
				temp = gsl_vector_ulong_get(pdata, count2);
				gsl_vector_ulong_set(pgraydata, count2,
						gsl_vector_ulong_get(pgrayindexes, temp));
//			gsl_vector_ulong_set(pgraydata,count2,gsl_vector_ulong_get(pdata,count2));   //find the gray
			}
			modulator(pdata, ptransmitted, psymbolconstellation); //ptransmitted stores the complex symbol

			channel_generator(pH, pr);
			sigman = sqrt((double)(1) / (2*SNR[count]*(float)(log2(M)))); /* corresponding noise standard deviation per dimension */
			noise_generator(pnoise, pr, sigman);
            gsl_blas_zcopy(pnoise, preceived);
			gsl_blas_zgemv(CblasNoTrans, one, pH, ptransmitted, one, preceived);

#ifdef RUNFCSD
			start = clock();
//			FCSD_CPU(preceived, pH, psymbolconstellation, SNR[count], symOut, durationKernel_CPU);

			end = clock();
	        Time_FCSD+=(double)(end-start);
#endif
//			cudaProfilerStart();
//			cudaDeviceReset();
#ifdef RUNCSVR_MIMO
			int start_M=0; //cold start
			int selection_M=0; //selection method WSS1D
			start = clock();
//			SVR_DETECTOR(preceived, pH, SNR[count], symOut, start_M,  selection_M);
			end = clock();
			Time_CSVR+=(double)(end-start);
#endif

//			demodulator_CPU(symOut, psymbolconstellation, pgrayindexes, poutput); //poutput is the graycode of the symbols
			iteration = iteration + 1;
			//check symbol error CPU version
#ifndef DEBUG
			for (count3 = 0; count3 < Nt; count3++) {
				if (fabs(gsl_vector_complex_get(symOut,count3).dat[0]-gsl_vector_complex_get(ptransmitted,count3).dat[0])>epsilon
						||fabs(gsl_vector_complex_get(symOut,count3).dat[1]-gsl_vector_complex_get(ptransmitted,count3).dat[1])>epsilon){
					symbolError=symbolError+1;
				}
			}
#endif
			//check symbol error GPU version
//			for (count3 = 0; count3 < Nt; count3++) {
//				if (symOut_cu[count3].x!=ptransmitted_cu[count3].x||symOut_cu[count3].y!=ptransmitted_cu[count3].y) {
//					symbolError=symbolError+1;
//				}
//			}
			//check bit error
#ifndef DEBUG
			for (count4 = 0; count4 < Nt; count4++) {
				errors = binaryerrors(gsl_vector_ulong_get(pgraydata, count4),
						gsl_vector_ulong_get(poutput, count4), M);
//        	 errors=1;
				BitError += errors;
			}
#endif
//			printf("the iteration time is %d\n", iteration);
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

//		} while ((symbolError < minSymbolError) || (iteration < miniteration));
		}while(iteration<1e4);
		time2=clock();
		Time_CSVR_aver[count]=Time_CSVR/(iteration*(CLOCKS_PER_SEC));//average time cost of CSVR
		Time_SVM_aver[count]=Time_SVM/(iteration*(CLOCKS_PER_SEC)); //average time cost of SVM
		Time_FCSD_aver[count]=Time_FCSD/(iteration*(CLOCKS_PER_SEC)); //average time cost of FCSD
		duration_total[count]=(double)(time2-time1)/(CLOCKS_PER_SEC); //total duration time of communication line
		printf(" SNR point %d finished\n", 2*count);
		printf("Iteration time at this SNR point is %f\n",(double)(iteration));
		BitErrorRate = (double)(BitError)/(double)((double)(iteration)*Nt*log2(M));
		SymbolErrorRate=(double)(symbolError)/(double)((double)(iteration)*Nt);
		if(count==begin){
			fprintf(cfile1, "BER: ");
			fprintf(cfile2, "SER: ");
		}
		fprintf(cfile1, "%0.20f ", BitErrorRate);
		fprintf(cfile2, "%0.20f ", SymbolErrorRate);
		if(count==SNRnum){
			fprintf(cfile1, "\n");
			fprintf(cfile2, "\n");
		}


	}

   fprintf(cfile3, "average time of CSVR is: ");
   for(count1=0;count1<SNRnum;count1++){
	   fprintf(cfile3, "%0.20f ", Time_CSVR_aver[count1]);
   }
   fprintf(cfile3, "\n");

   fprintf(cfile3, "average time of SVM is: ");
   for(count1=0;count1<SNRnum;count1++){
	   fprintf(cfile3, "%0.20f ", Time_SVM_aver[count1]);
   }
   fprintf(cfile3, "\n");

   fprintf(cfile3, "average time of FCSD is: ");
   for(count1=0;count1<SNRnum;count1++){
	   fprintf(cfile3, "%0.20f ", Time_FCSD_aver[count1]);
   }
   fprintf(cfile3, "\n");

   fprintf(cfile3, "total duration of communication line is: ");
   for(count1=0;count1<SNRnum;count1++){
	   fprintf(cfile3, "%0.20f ", duration_total[count1]);
   }
   fprintf(cfile3, "\n");

   fprintf(cfile4, "average computational cost of CSVR: ");
   for(count1=0;count1<SNRnum; count1++){

   }
   fprintf(cfile4, "\n");

   fprintf(cfile4,"average iteration time of CSVR: ");
   for(count1=0;count1<SNRnum;count1++){

   }
   fprintf(cfile4, "\n");

#ifdef DEBUG
   fprintf(cfile2, "\n");
   fprintf(cfile2, "The program terminated normally\n");
#endif

	gsl_rng_free(pr);
	pr = NULL;
	gsl_vector_ulong_free(pgraycode);
	pgraycode = NULL;
	gsl_vector_ulong_free(pgrayindexes);
	pgrayindexes = NULL;
	gsl_vector_complex_free(psymbolconstellation);
	psymbolconstellation = NULL;
	free(SNR);
	SNR = NULL;
	printf("the whole program finishes successfully\n");
	fclose(cfile1);
	cfile1=NULL;
	fclose(cfile2);
	cfile2=NULL;
	fclose(cfile3);
	cfile3=NULL;
	fclose(cfile4);
	cfile4=NULL;
	return (0);

}
