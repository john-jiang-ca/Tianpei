/*
 ============================================================================
 Name        : CommunicationLink.c
 Author      : Tianpei Chen
 Version     : beta 1.0
 Copyright   : Your copyright notice
 Description : Communication Link
  This platform is the general software testbed for the performance test of
LS-MIMO detectors proposed. The simulation results are presented in the form
of bit error rate (BER), symbol error rate (SER) and frame error rate (FER)
versus signal to noise ratio (SNR) in dB.

  The SNR considered here is the average receive SNR (general) or symbol SNR (Dejelili"
Thesis). At each SNR point, Monte-Carlo (MC) simulation is performed in the step
of channel realizations. The MC simulation stops until a certain number of channel
realization is achieved as well as a certain number of symbol errors are accumulated.

  The simulation results are recorded in the output files.
 ============================================================================
 */


#include "commonSettings.h"

void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, int M); //generate transmit index
void grayencoder (gsl_vector_ulong *pgraydata); // gray encoder
void symbolconstellation(gsl_vector_complex *psymbolconstellation, double pav); //symbol constellation generator
void modulator (gsl_vector_ulong *pdata,  gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,
		gsl_vector_complex *ptransmitted, gsl_vector_ulong *pgrayInput); //modulation
void error_channel_generator (gsl_matrix_complex *pH, gsl_rng *pr, double sigmaerr); //channel estimation error generator
void corr_matrix_generator (gsl_matrix_complex *pRt, gsl_matrix_complex *pRr, double paraT, double paraR);  //correlation matrix generator
void channel_generator (gsl_matrix_complex *pH, gsl_rng *pr); //random channel generator
void noise_generator (gsl_vector_complex *pnoise, gsl_rng *pr, double noiseV); //AWGN generator
void MMSE(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr, double pav, int M,  gsl_vector_complex *psymOut); //MMSE detector
void symErrorCheck(gsl_vector_complex *ptransmit, gsl_vector_complex *psymOut, int *ErrorIndex, int *symError_sub, int *frameError_sub);
//find the error symbol"s index and calculate the number of symbol error
void demodulator(gsl_vector_complex *psymOut, gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,
		int *ErrorIndex,  gsl_vector_ulong *pgrayOut);
// demodulation for the error symbol
void binaryerrors (gsl_vector_ulong *pgrayOut, int *ErrorIndex, gsl_vector_ulong *pgrayInput,  int M, int *bitError_sub);




int main(void) {
	int Nr=receiveAntennas;
	int Nt=transmitAntennas;
	int N=floor(sqrt(Nr+(1/4)*pow((Nr-Nt),2.0))-(1/2)*(Nr-Nt));  //the number of antennas that are chosen in the channel partition stage
	int M=symConstellationSize;
    int SNR_tmp;
    double pav=(double)1/((double)Nt); //the average power of transmit symbol
    clock_t start, end;
    FILE *pfile;
    pfile=fopen(fileName, "a");
    fprintf(pfile, "==============================================================================\n");
    fprintf(pfile, "This the output file for software testbed in Large-Scale MIMO (LS-MIMO) system\n");
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "SYSTEM CONFIGURATION\n");
    fprintf(pfile, "this is for %d times %d MIMO with %d QAM modulation\n", Nr, Nt, M);
    fprintf(pfile, "the number of antennas chosen in the channel partition is %d\n", N);
    if (Corr_Ind==1){
    	fprintf(pfile, "consider the channel correlation\n");
    	fprintf(pfile, "correlation parameter Rr=%g, Rt=%g\n",  Rr,  Rt);
    }else{
    	fprintf(pfile, "do not consider spatial correlation\n");
    }
    if(Est_Ind==1){
    	fprintf(pfile, "consider channel estimation error\n");
    	fprintf(pfile, "Estimation error parameter is %g",  gammasq);
    }else{
    	fprintf(pfile, "do not consider the channel estimation error\n");
    }
    fprintf(pfile, "the average transmit symbol power is %g\n", pav);
    SNR_tmp=Start_SNR;
    fprintf(pfile, "SNR (dB) are\n");
    while(SNR_tmp<=End_SNR){
    	fprintf(pfile, "%d, ", SNR_tmp );
    	SNR_tmp+=Step_SNR;
    }
    fprintf(pfile, "\n");
    fprintf(pfile, "the minimum channel realization is %g\n", (double) minChannelRealizations);
    fprintf(pfile, "the minimum symbol errors accumulated is %g\n", (double) minSymErrors);
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "\n");
    fprintf(pfile, "\n");

    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "SYSTEM OUTPUT\n");
    fprintf(pfile, "SNR \t Channel Realization \t frame errors \t symbol errors \t bit errors \t FER \t SER \t BER \t Operation time\n");
    fclose(pfile);

    gsl_vector_ulong *pgraydata = gsl_vector_ulong_calloc (M); //gray code book
    gsl_vector_complex *psymbolconstellation = gsl_vector_complex_calloc (M); //symbol constellation
    gsl_vector_ulong *pdata = gsl_vector_ulong_calloc (Nt); //the indexes of the transmit symbol vector
    gsl_vector_ulong *pgrayInput = gsl_vector_ulong_calloc (Nt);// the gray code of the transmit data
    gsl_vector_complex *ptransmitted = gsl_vector_complex_calloc (Nt); //the transmitted symbol vector
    gsl_matrix_complex *pH = gsl_matrix_complex_calloc (Nr, Nt);   //the channel matrix
    gsl_matrix_complex *pH_tmp=gsl_matrix_complex_calloc(Nr, Nt);
    gsl_matrix_complex *pHest = gsl_matrix_complex_calloc (Nr, Nt);   //estimation error matrix
    gsl_matrix_complex *pRr = gsl_matrix_complex_calloc (Nr, Nr);  //receive spatial correlation matrix
    gsl_matrix_complex *pRt = gsl_matrix_complex_calloc (Nt, Nt);  //transmit spatial correlation matrix
    gsl_vector_complex *pnoise = gsl_vector_complex_calloc (Nr);  //AWGN vector
    gsl_vector_complex *preceived = gsl_vector_complex_calloc (Nr);  //the received symbol vector
    gsl_vector_complex *psymOut=gsl_vector_complex_calloc(Nt); //the detected of the symbol vector
    int *ErrorIndex_V=(int*)malloc(sizeof(int)*Nt);    //the indexes of the erroneous symbol
    gsl_vector_ulong *pgrayOut;   //output gray code for erroneous symbols
    gsl_rng *pr;    //random number generator
	int seed = time (NULL);
	const gsl_rng_type *pT;
	pT = gsl_rng_default;
	pr = gsl_rng_alloc (pT);
    gsl_rng_set (pr,seed);
    gsl_complex alpha, beta;
    GSL_SET_COMPLEX(&alpha, 1,0);
    GSL_SET_COMPLEX(&beta, 0, 0);
    int frameError, symError, bitError;  //accumulated frame errors, symbol errors and bit errors in one SNR point
    int *frameError_sub, *symError_sub, *bitError_sub; //frame errors, symbol errors and bit errors in one channel realization
    frameError_sub=(int*)malloc(sizeof(int));
    symError_sub=(int*)malloc(sizeof(int));
    bitError_sub=(int*)malloc(sizeof(int));
    double  noiseV, snr,sigmaerr;
    double FER, SER, BER;   //frame error rate, symbol error rate and bit error rate
    int Realizations=0;   //the number of channel realizations
    SNR_tmp=Start_SNR;
    grayencoder (pgraydata);   //generate gray code book
    int count;
#ifdef DEBUG
    printf("the gray data are\n");
    for (count=0;count<M;count++){
    	printf("%d, ", gsl_vector_ulong_get(pgraydata, count));
    }
    printf("\n");
#endif
    symbolconstellation(psymbolconstellation, pav);  //generate the symbol constellation
//    gsl_complex constemp;
//    for (count=0;count<M;count++){
//    	double a=gsl_vector_complex_get(psymbolconstellation, count).dat[0];
//    	double b=gsl_vector_complex_get(psymbolconstellation, count).dat[1];
//    }
    if(Corr_Ind==1){
    	corr_matrix_generator (pRt, pRr,  Rt,  Rr);  //generate spatial correlation matrix
    }
    while (SNR_tmp<=End_SNR){
    	frameError=0;
    	symError=0;
    	bitError=0;
    	FER=0;
    	SER=0;
    	BER=0;
    	Realizations=0;
    	pfile=fopen(fileName, "a");
    	start=clock();
    	while (symError<minSymErrors||Realizations<minChannelRealizations){
    		data_generator(pdata, pr, M); //generate the random index of the data to be transmitted
#ifdef DEBUG
    		printf("Test of pdata\n");

    		for ( count=0;count<Nt; count++){
    			printf("%d, ", gsl_vector_ulong_get(pdata, count));
    		}
    		printf("\n");
#endif
    		unsigned long temp;
    		for(count=0;count<Nt;count++){
    			temp=gsl_vector_ulong_get(pdata, count);
    		 gsl_vector_ulong_set(pgrayInput, count, gsl_vector_ulong_get(pgraydata, temp));
//    		 gsl_vector_complex_set(ptransmitted, count, gsl_vector_complex_get(psymbolconstellation, temp));
    		}
//    		modulator (pdata,  psymbolconstellation, pgraydata, ptransmitted, pgrayInput); //generate the corresponding transmit
    		//symbol vector and gray code vector
    		channel_generator (pH, pr);
    		if(Corr_Ind==1){
    			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, pRr, pH, beta, pH_tmp);   //add spatial correlation (receive side)
    			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, pH_tmp, pRt, beta, pH);   //add spatial correlation (transmit side)
    		}
    		snr=pow(10,(SNR_tmp/10)); //SNR in decimal
    		noiseV=1/snr;     //noise variance
    		if(Est_Ind==1){
    			sigmaerr = sqrt (gammasq*pav/(2.0*snr));
    			error_channel_generator (pH_tmp, pr, sigmaerr);    //add channel estimation error
    			gsl_matrix_complex_memcpy(pHest, pH);
    			gsl_matrix_complex_add(pHest,pH_tmp);   //generate imperfect channel estimation
    		}
    		noise_generator (pnoise, pr, noiseV);   //generate noise vector
    		gsl_vector_complex_memcpy(preceived, pnoise);
    		gsl_blas_zgemv(CblasNoTrans, alpha, pH, ptransmitted, alpha, preceived);  //generate  receive signal vector
    		if(Est_Ind==1){
    		MMSE(preceived, pHest, snr/(double)(Nt), pav,  M,  psymOut); //detection with imperfect CSI
    		}else{
    	    MMSE(preceived, pH, snr/(double)(Nt), pav,  M,  psymOut); //detection with perfect CSI
    		}
    		//MMSE detector
    		symErrorCheck(ptransmitted, psymOut, ErrorIndex_V, symError_sub, frameError_sub);   //check symbol error
    		Realizations++;
    		if (*symError_sub==0){
    			continue;
    		}
    		int *ErrorIndex=(int*)malloc(sizeof(int)*(*symError_sub));
    		for (count=0;count<(*symError_sub); count++){
    			ErrorIndex[count]=ErrorIndex_V[count];
    		}
    		pgrayOut=gsl_vector_ulong_calloc(symError_sub[0]);
    		demodulator(psymOut, psymbolconstellation, pgraydata, ErrorIndex,  pgrayOut);
    		binaryerrors (pgrayOut, ErrorIndex, pgrayInput,  M,  bitError_sub);
    		symError+=*symError_sub;
    		frameError+=*frameError_sub;
    		bitError+=*bitError_sub;
    		free(ErrorIndex);

    	}
        end=clock();
        FER=(double)frameError/((double)(Realizations));  //calculate frame error rate
        SER=(double)symError/((double)(Realizations*Nt)); //calculate symbol error rate
        BER=(double)bitError/((double)(Realizations*Nt*ceil(log2(M)))); //calculate bit error rate
    	printf("SNR=%d, Realization=%d, FER=%g, SER=%g, BER=%g, OperationTime=%g s\n", SNR_tmp, Realizations, FER, SER, BER,
    			(end-start)/(double)CLOCKS_PER_SEC);
    	fprintf(pfile, "%d \t %d \t %d \t %d \t %d \t %g \t %g \t %g \t %g\n", SNR_tmp, Realizations, frameError, symError, bitError, FER, SER, BER, (end-start)/(double)CLOCKS_PER_SEC);
    	fclose(pfile);
    	SNR_tmp+=Step_SNR;

    }
    pfile=fopen(fileName, "a");
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "The whole program ends successfully!!\n");
    fprintf(pfile, "==============================================================================\n");
    fclose(pfile);
    printf("the whole program ends successfully!!\n");
    free(pfile);
    gsl_vector_ulong_free(pgraydata);
    gsl_vector_complex_free(psymbolconstellation);
    gsl_vector_ulong_free(pdata);
    gsl_vector_ulong_free(pgrayInput);
    gsl_vector_complex_free(ptransmitted);
    gsl_matrix_complex_free(pH);
    gsl_matrix_complex_free(pH_tmp);
    gsl_matrix_complex_free(pHest);
    gsl_matrix_complex_free(pRr);
    gsl_matrix_complex_free(pRt);
    gsl_vector_complex_free(pnoise);
    gsl_vector_complex_free(preceived);
    gsl_vector_complex_free(psymOut);
    gsl_vector_ulong_free(pgrayOut);
    gsl_rng_free(pr);
    free(frameError_sub);
    free(symError_sub);
    free(bitError_sub);
    free(ErrorIndex_V);
	return 0;
}
