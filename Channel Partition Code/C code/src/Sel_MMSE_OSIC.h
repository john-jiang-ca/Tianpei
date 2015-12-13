/*
 * Sel_MMSE_OSIC.h
 * Notice  : the input signal to noise ratio is average receive
 * signal to noise ratio
 *  Created on: Dec 12, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "diversity_max_selection.h"
#include "MMSE_OSIC.h"
#include "linkedList.h"
#include "fullfact.h"
#ifndef SEL_MMSE_OSIC_H_
#define SEL_MMSE_OSIC_H_
void Sel_MMSE_OSIC(gsl_vector_complex *preceive, gsl_matrix_complex *pH,
		double SNR, double pav, int N, gsl_vector_complex *psymbolconstellation,
		gsl_vector_complex *psymOut){
	int Nr=pH->size1;
	int Nt=pH->size2;
	int M=psymbolconstellation->size;
	int listSize=pow(M, N);   //the number of vector candidates in the list
	gsl_matrix_complex *pH1, *pH2;
    gsl_vector_complex *preceiveSub, *preceiveSubtmp, *psymOut_sub1, *psymOut_sub2,
    *psymOut_sub2tmp;
    gsl_matrix_complex *psymOut_sub2_M;
    gsl_matrix *List;
    gsl_vector *subList;
    gsl_vector *EuclideanV;
    struct Node *index1_head, *index2_head;
	pH1=gsl_matrix_complex_calloc(Nr, N);
	pH2=gsl_matrix_complex_calloc(Nr, Nt-N);
    preceiveSub=gsl_vector_complex_calloc(Nr);
    preceiveSubtmp=gsl_vector_complex_calloc(Nr);
    psymOut_sub1=gsl_vector_complex_calloc(N);
    psymOut_sub2=gsl_vector_complex_calloc(Nt-N);
    psymOut_sub2tmp=gsl_vector_complex_calloc(Nt-N);
    psymOut_sub2_M=gsl_matrix_complex_calloc(Nt-N, listSize);
    subList=gsl_vector_calloc(N);
    List=gsl_matrix_calloc(N,listSize);
    EuclideanV=gsl_vector_calloc(listSize);
    double Euclidean;
    int count, count1;
    gsl_complex alpha, beta1, beta2;
    GSL_SET_COMPLEX(&alpha, 1, 0);
    GSL_SET_COMPLEX(&beta1, 0, 0);
    GSL_SET_COMPLEX(&beta2, -1, 0);
    double snr=SNR/(double)Nt;   //symbol signal to noise ratio
    int k;
    int index;
    fullfact(N, M, List);   //generate all the possible sub symbol vector candidate x1
    diversity_max_selection(pH, N, snr, index1_head, index2_head, pH1, pH2);
    for (count=0;count<listSize; count++){  //perform MMSE-OSIC detection to d2
    	gsl_matrix_get_col(subList, List, count);
    	for (count1=0;count1<N;count1++){  //construct sub symbol vector
    		gsl_vector_complex_set(psymOut_sub1, count1,
    				gsl_vector_complex_get(psymbolconstellation,
    						gsl_vector_get(subList, count1)));
    	}
    	gsl_blas_zgemv(CblasNoTrans, alpha, pH1, psymOut_sub1, beta1, preceiveSubtmp);
    	gsl_vector_complex_memcpy(preceiveSub, preceive);
    	gsl_vector_complex_sub(preceiveSub, preceiveSubtmp);
    	MMSE_OSIC(preceiveSub, pH2, snr, pav, M, psymOut_sub2tmp);
    	gsl_matrix_complex_set_col(psymOut_sub2_M, count, psymOut_sub2tmp);
    	gsl_blas_zgemv(CblasNoTrans, alpha, pH2, psymOut_sub2tmp, beta2, preceiveSub);
    	Euclidean=gsl_blas_dznrm2(preceiveSub);
    	gsl_vector_set(EuclideanV, count, Euclidean);

    }
    //reconstruct psymOut
    k=gsl_vector_min_index(EuclideanV);
   	gsl_matrix_get_col(subList, List, k);
    	for (count1=0;count1<N;count1++){  //construct sub symbol vector
    		gsl_vector_complex_set(psymOut_sub1, count1,
    				gsl_vector_complex_get(psymbolconstellation,
    						gsl_vector_get(subList, count1)));
    	}

    	gsl_matrix_complex_get_col(psymOut_sub2, psymOut_sub2_M, k);

    	for (count1=0;count1<N;count1++){
         index=get(0, index1_head);
         index1_head=index1_head->next;
         gsl_vector_complex_set(psymOut, index,
        		 gsl_vector_complex_get(psymOut_sub1,count1));
    	}
    	for (count1=0;count1<(Nt-N);count1++){
    		index=get(0, index2_head);
    		index2_head=index2_head->next;
            gsl_vector_complex_set(psymOut, index,
           		 gsl_vector_complex_get(psymOut_sub2,count1));
    	}

	gsl_matrix_complex_free(pH1);
	gsl_matrix_complex_free(pH2);
	gsl_vector_complex_free(preceiveSub);
	gsl_vector_complex_free(preceiveSubtmp);
	gsl_vector_complex_free(psymOut_sub1);
	gsl_vector_complex_free(psymOut_sub2);
	gsl_vector_complex_free(psymOut_sub2tmp);
	gsl_matrix_complex_free(psymOut_sub2_M);
	gsl_matrix_free(List);
	gsl_vector_free(subList);
	gsl_vector_free(EuclideanV);

	return;
}




#endif /* SEL_MMSE_OSIC_H_ */
