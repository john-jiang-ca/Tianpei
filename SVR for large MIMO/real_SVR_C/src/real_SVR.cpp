/*
 ============================================================================
 Name        : real_SVR.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include"Public.h"
#include"WSS2D_1Dsolver.h"
#include"MMSE_detec.h"
int main(void) {
	gsl_matrix *pH=gsl_matrix_calloc(Nr,Nt);
	gsl_vector *symTransmitted=gsl_vector_calloc(Nt);
	gsl_vector *symReceived=gsl_vector_calloc(Nr);
	gsl_vector *noise=gsl_vector_calloc(Nr);
    gsl_vector *SNRd=gsl_vector_calloc(SNRnum);
    gsl_vector *noiseVariance=gsl_vector_calloc(SNRnum);
	gsl_vector *symOut=gsl_vector_calloc(Nt);
	gsl_vector *symbolconstellation=gsl_vector_calloc(M);

//	gsl_vector *sigma1=gsl_vector_calloc(Nr);

    int TransmitNum[Nt];
    int SNR[SNRnum];
    double symError;
    double symErrorRate[SNRnum];
    int channel_realization;
    double d;

//*iteration=0;
	unsigned long seed=0L;
	const gsl_rng_type *pT;
	pT=gsl_rng_default;
    gsl_rng *pr=gsl_rng_alloc(pT);
    seed=time(NULL);
    gsl_rng_set(pr,seed);
//    double SNR_current;
printf("the program begin!!");

    int count, count1, count2;
    printf("the SNR point are: \n");
    for(count=0;count<SNRnum;count++){
    	SNR[count]=2*count;
//    	printf("%d ", SNR[count]);
    	gsl_vector_set(SNRd, count, pow(10, (double)SNR[count]/(double)10));
    	gsl_vector_set(noiseVariance,count, (double)1/gsl_vector_get(SNRd,count));
    	printf("%d %f ", SNR[count], gsl_vector_get(noiseVariance,count));
    }
printf("\n");
//define the symbol constellation
d=sqrt((double) 3/(Nt*(pow(M,2)-1)));
for(count1=0;count1<M;count1++){
	gsl_vector_set(symbolconstellation,count1, -(M-1)*d+count1*2*d);
	printf("the %d th symbol is %f\n ", count1+1, gsl_vector_get(symbolconstellation,count1));
}
printf("\n");
//int begin=11;//define start SNR point
int aver_iter=0;
//initialize the file to store data
    FILE *cfile1, *cfile2, *cfile3;
    cfile1=fopen(Iteration, "a");
    cfile2=fopen(SER, "a");
    cfile3=fopen(BER, "a");

    //file1
    fprintf(cfile1, "\n");
#ifdef TEST
    fprintf(cfile1, "test\n");
#endif

#ifndef TEST
    fprintf(cfile1, "formal\n");
#endif
#ifdef SVR
    fprintf(cfile1, "SVR\n");
#endif
#ifdef MMSE
    fprintf(cfile1, "MMSE\n");
#endif
    fprintf(cfile1, "hyperparameter settings: \n");
    fprintf(cfile1, "epsilon %0.20f, tol %0.20f, C %f\n", epsilon, tol, C);
    fprintf(cfile1, "%d X %d MIMO system  %d PAM modulation\n", Nr, Nt, M);
//    fprintf(cfile1, "SNR from %d to %d\n", 2*begin, 2*SNRnum);
    fprintf(cfile1, "The SNR point is: %d %f to %d %f\n", SNR[begin],
    		gsl_vector_get(SNRd, begin), SNR[SNRnum-1],gsl_vector_get(SNRd, SNRnum-1) );
//    for(count1=0;count1<SNRnum-begin+1;count1++){
//    	fprintf(cfile1,"%f ", gsl_vector_get(SNRd,count1));
//    }
    fprintf(cfile1,"\n");
    fprintf(cfile1, "The average iteration time are: \n");

     //file 2
    fprintf(cfile2, "\n");
#ifdef TEST
    fprintf(cfile2, "test\n");
#endif

#ifndef TEST
    fprintf(cfile2, "formal\n");
#endif
#ifdef SVR
    fprintf(cfile2, "SVR\n");
#endif
#ifdef MMSE
    fprintf(cfile2, "MMSE\n");
#endif
    fprintf(cfile2, "hyperparameters settings: \n");
    fprintf(cfile2, "epsilon %0.20f, tol %0.20f, C %f\n", epsilon, tol, C);
    fprintf(cfile2, "%d X %d MIMO system  %d PAM modulation\n", Nr, Nt, M);
//    fprintf(cfile2, "SNR from %d to %d\n", 2*begin, 2*SNRnum);
    fprintf(cfile2, "The SNR point is: %d %f to %d %f\n", SNR[begin],
    		gsl_vector_get(SNRd, begin), SNR[SNRnum-1],gsl_vector_get(SNRd, SNRnum-1) );
//    for(count1=0;count1<SNRnum-begin+1;count1++){
//    	fprintf(cfile2,"%f ", gsl_vector_get(SNRd,count1));
//    }
    fprintf(cfile2,"\n");
    fprintf(cfile2, "symbol error rate are: \n");
    //file 3

    fprintf(cfile3, "\n");
#ifdef TEST
    fprintf(cfile3, "test\n");
#endif

#ifndef TEST
    fprintf(cfile3, "formal\n");
#endif
#ifdef SVR
    fprintf(cfile3, "SVR\n");
#endif
#ifdef MMSE
    fprintf(cfile3, "MMSE\n");
#endif
    fprintf(cfile3, "hyperparameter settings: \n");
    fprintf(cfile3, "epsilon %0.20f, tol %0.20f, C %f\n", epsilon, tol, C);
    fprintf(cfile3, "\n");
    fprintf(cfile3, "%d X %d MIMO system  %d PAM modulation\n", Nr, Nt, M);
//    fprintf(cfile3, "SNR from %d to %d\n", 2*begin, 2*SNRnum);
    fprintf(cfile3, "The SNR point is: %d %f to %d %f\n", SNR[begin],
    		gsl_vector_get(SNRd, begin), SNR[SNRnum-1],gsl_vector_get(SNRd, SNRnum-1) );
//    for(count1=0;count1<SNRnum-begin+1;count1++){
//    	fprintf(cfile3,"%f ", gsl_vector_get(SNRd,count1));
//    }
    fprintf(cfile3,"\n");
    fprintf(cfile3, "Bit error rate are: \n");
    fclose(cfile1);
    fclose(cfile2);
    fclose(cfile3);

for(count=begin;count<SNRnum;count++) {
symError=0;
channel_realization=0;
aver_iter=0;
cfile1=fopen(Iteration, "a");
cfile2=fopen(SER, "a");
cfile3=fopen(BER, "a");
while(channel_realization<1e5||symError<300){
	gsl_vector *lamida=gsl_vector_calloc(Nr);
	double *G=(double*)malloc(sizeof(double));
	double *Theta=(double*)malloc(sizeof(double));
	double *MSE=(double*)malloc(sizeof(double));
    int *iteration=(int*)malloc(sizeof(int));
    //Generate random channel propagation matrix and AWGN noise

    for(count1=0;count1<Nr;count1++){
    	for(count2=0;count2<Nt;count2++){
    		gsl_matrix_set(pH, count1,count2, gsl_ran_gaussian(pr, 1));
    	}
    	gsl_vector_set(noise, count1, gsl_ran_gaussian(pr, sqrt(gsl_vector_get(noiseVariance,count))));
    }



    //modulation
#ifdef DEBUG
printf("the random data generated is: \n");
#endif
    for(count1=0;count1<Nt;count1++){
        TransmitNum[count1]=gsl_rng_uniform_int(pr, M);
    	gsl_vector_set(symTransmitted, count1, gsl_vector_get(symbolconstellation,TransmitNum[count1]));
#ifdef DEBUG
    	printf("%d ", TransmitNum[count1]);
#endif
    }
#ifdef DEBUG
    printf("\n");
#endif

    //AWGN channel modeling
    gsl_blas_dcopy(noise, symReceived);
     gsl_blas_dgemv(CblasNoTrans, 1, pH, symTransmitted, 1, symReceived);
#ifdef DEBUG
     printf("the level of received symbol vector\n");
     for(count2=0;count2<Nr;count2++){
    	 printf("%f ", gsl_vector_get(symReceived,count2));
     }
     printf("\n");
#endif
     //detector
#ifdef SVR
     WSS2D_1Dsolver(pH, symReceived, gsl_vector_get(SNRd,count), symbolconstellation,
    		 symOut, lamida, Theta, G, iteration, MSE);
#endif


#ifdef MMSE
     MMSE_detec(pH, symReceived, gsl_vector_get(SNRd,count), symbolconstellation,symOut,  MSE);
#endif
#ifdef DEBUG
     printf("the  output symbol is\n");
     for(count1=0;count1<Nt;count1++){
    	 printf("%f ", gsl_vector_get(symOut,count1));
     }
     printf("the  input symbol is\n");
     for(count1=0;count1<Nt;count1++){
    	 printf("%f ", gsl_vector_get(symTransmitted,count1));
     }
#endif
     //calculate symbol errors
     for(count1=0;count1<Nt;count1++){
    	 if(fabs(gsl_vector_get(symOut, count1)-gsl_vector_get(symTransmitted,count1))>1e-5)
    	 {
    		 symError++;
    	 }
     }
     aver_iter+=*iteration;
     channel_realization++;
     free(G);
     free(Theta);
     free(MSE);
     free(iteration);
     gsl_vector_free(lamida);
}
symErrorRate[count]=(double)symError/(double)(Nt*channel_realization);
aver_iter=aver_iter/channel_realization;
fprintf(cfile1, "%d ", aver_iter);
fprintf(cfile2, "%0.20f ", symErrorRate[count]);
//#ifdef DEBUG
printf("SNR point %d end", SNR[count]);
printf("\n");
//#endif
fclose(cfile1);
fclose(cfile2);
fclose(cfile3);
}
#ifdef DEBUG
printf("the symbol error rate is \n");
for (count=0;count<SNRnum;count++){

	printf("%f ", symErrorRate[count]);
}
printf("\n");
#endif

gsl_matrix_free (pH);
gsl_vector_free (symTransmitted);
gsl_vector_free (symReceived);
gsl_vector_free (noise);
gsl_vector_free(SNRd);
gsl_vector_free (noiseVariance);
gsl_vector_free (symOut);
gsl_vector_free(symbolconstellation);

//gsl_vector_free(sigma1);
gsl_rng_free(pr);
cfile1=fopen(Iteration,"a");
cfile2=fopen(SER, "a");
cfile3=fopen(BER, "a");
fprintf(cfile1, "\n");
fprintf(cfile2, "\n");
fprintf(cfile3, "\n");
fprintf(cfile1, "The program end");
fprintf(cfile2, "The program end");
fprintf(cfile3, "The program end");
fclose(cfile1);
fclose(cfile2);
fclose(cfile3);
printf("The whole program end!!");
return 0;


}
