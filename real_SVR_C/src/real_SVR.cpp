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
//#include"WSS1D_2Dsolver.h"

int main(void) {
	gsl_matrix *pH=gsl_matrix_calloc(Nr,Nt);
//	gsl_matrix *kernel=gsl_matrix_calloc(Nt, Nt);
	gsl_vector *symTransmitted=gsl_vector_calloc(Nt);
	gsl_vector *symReceived=gsl_vector_calloc(Nr);
	gsl_vector *noise=gsl_vector_calloc(Nr);
    gsl_vector *SNRd=gsl_vector_calloc(SNRnum);
    gsl_vector *noiseVariance=gsl_vector_calloc(SNRnum);
	gsl_vector *symOut=gsl_vector_calloc(Nt);
	gsl_vector *symbolconstellation=gsl_vector_calloc(M);
	gsl_vector *lamida=gsl_vector_calloc(Nr);
	gsl_vector *sigma=gsl_vector_calloc(Nr);
	double G[Size];
	double Theta[Size];
	double Theta_tmp;
	double lamida_tmp;
    int TransmitNum[Nt];
    int SNR[SNRnum];
    int First, Second;
    double symError;
    double symErrorRate;
    int iteration;
    int channel_realization;
    double d;
//    double SNR_current;
printf("the program begin!!");

    int count, count1, count2;
    printf("the SNR point are: \n");
    for(count=0;count<SNRnum;count++){
    	SNR[count]=2*count;
    	printf("%d ", SNR[count]);
    	gsl_vector_set(SNRd, count, pow(10, (double)SNR[count]/10));
    }

    for(count=0;count<SNRnum;count++){
    	gsl_vector_set(noiseVariance,count, (double)1/gsl_vector_get(SNRd,count));
    }

    FILE *cfile1, *cfile2, *cfile3;
    cfile1=fopen(Iteration, "a");
    cfile2=fopen(SER, "a");
    cfile3=fopen(BER, "a");

for(count=0;count<SNRnum;count++) {

    //Generate random number
	const gsl_rng_type *pT;
	pT=gsl_rng_default;
    gsl_rng *pr=gsl_rng_alloc(pT);
    for(count1=0;count1<Nr;count1++){
    	for(count2=0;count2<Nt;count2++){
    		gsl_matrix_set(pH, count1,count2, gsl_ran_gaussian(pr, 1));
    	}
    	gsl_vector_set(noise, count1, gsl_ran_gaussian(pr, gsl_vector_get(noiseVariance,count)));
    }



    //modulation
    d=sqrt((double) 3/(Nt*(pow(M,2)-1)));
    for(count1=0;count1<M;count1++){
    	gsl_vector_set(symbolconstellation,count1, -(M-1)*d+count1*2*d);
    }
    for(count1=0;count1<Nt;count1++){
    	gsl_vector_set(symTransmitted, count1, TransmitNum[count1]);
    }

    //AWGN channel modeling
    gsl_blas_dcopy(noise, symReceived);
     gsl_blas_dgemv(CblasNoTrans, 1, pH, symTransmitted, 1, symReceived);
     //detector
#ifndef DEBUG
     WSS2D_1Dsolver();
#endif
     //calculate symbol errors
     for(count1=0;count1<Nt;count1++){
    	 if(fabs(gsl_vector_get(symOut, count1)-gsl_vector_get(symTransmitted,count1))>1e-5)
    	 {
    		 symError+=1;
    	 }
     }

}

fclose(cfile1);
fclose(cfile2);
fclose(cfile3);
return 0;


}
