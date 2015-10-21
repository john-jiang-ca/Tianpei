/*
 * MMSE_detec.h
 *
 *  Created on: Sep 2, 2015
 *      Author: Preston Chen
 */

#ifndef MMSE_DETEC_H_
#define MMSE_DETEC_H_

#include"Public.h"
void MMSE_detec(
		gsl_matrix *pH,
		gsl_vector *symReceived,
		double SNRd,
		gsl_vector *symConstellation,
		gsl_vector *symOut,
		double *MSE
		){
	//MMSE equalization
	gsl_matrix *G=gsl_matrix_calloc(Nt,Nt);
	gsl_permutation *p=gsl_permutation_calloc(Nt);
	gsl_matrix *LU=gsl_matrix_calloc(Nt,Nt);
    gsl_vector *y_tmp=gsl_vector_calloc(Nr);
    gsl_vector *symOut_tmp=gsl_vector_calloc(Nr);
	gsl_matrix *inverse=gsl_matrix_calloc(Nt,Nt);
	gsl_matrix *E=gsl_matrix_calloc(Nt, Nr);
    int *signum=(int*)malloc(sizeof(int));
	int count,count1,count2;
    gsl_matrix_set_identity(G);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, pH, pH, (double)1/(SNRd), G);
	gsl_matrix_memcpy(LU, G);
	gsl_linalg_LU_decomp(LU, p, signum);
	gsl_linalg_LU_invert(LU, p, inverse);
#ifdef DEBUG
	gsl_matrix *I=gsl_matrix_calloc(Nt,Nt);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverse, G, 0,I);
	printf("the identity matrix is\n");
	for(count1=0;count1<Nt;count1++){
		for(count2=0;count2<Nt;count2++){
			printf("%f ", gsl_matrix_get(I,count1,count2));
		}
		printf("\n");
	}
	gsl_matrix_free(I);
#endif

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, inverse, pH, 0, E);
	gsl_blas_dgemv(CblasNoTrans, 1, E, symReceived, 0, symOut);
#ifdef DEBUG
	printf("the unrounded output symbols are\n");
	for(count=0;count<Nt;count++){
		printf("%f ",gsl_vector_get(symOut,count));
	}
	printf("\n");
#endif

	//rounding
	double d=sqrt((double)3/(Nt*(pow(M,2)-1)));
	for(count1=0;count1<Nt;count1++){
		if(M==4){
			if(gsl_vector_get(symOut, count1)<=-2*d){
	         gsl_vector_set(symOut,count1, -3*d);
			}else if (gsl_vector_get(symOut ,count1)>-2*d&&gsl_vector_get(symOut, count1)<=0){
	 gsl_vector_set(symOut, count1, -d);
			}else if(gsl_vector_get(symOut,count1)>0&&gsl_vector_get(symOut, count1)<=2*d){
	gsl_vector_set(symOut, count1, d);
			}else if(gsl_vector_get(symOut ,count1)>2*d){
	gsl_vector_set(symOut, count1 , 3*d);
			}

		}else if(M==8){
			if(gsl_vector_get(symOut, count1)<=-6*d){
	         gsl_vector_set(symOut,count1, -7*d);
			}else if (gsl_vector_get(symOut ,count1)>-6*d&&gsl_vector_get(symOut, count1)<=-4*d){
	 gsl_vector_set(symOut, count1, -5*d);
			}else if(gsl_vector_get(symOut,count1)>-4*d&&gsl_vector_get(symOut, count1)<=-2*d){
	gsl_vector_set(symOut, count1, -3*d);
			}else if(gsl_vector_get(symOut ,count1)>-2*d&&gsl_vector_get(symOut, count1)<=0){
	gsl_vector_set(symOut, count1 , -d);
			}else if(gsl_vector_get(symOut, count1)>6*d){
		         gsl_vector_set(symOut,count1, 7*d);
				}else if (gsl_vector_get(symOut ,count1)>4*d&&gsl_vector_get(symOut, count1)<=6*d){
		 gsl_vector_set(symOut, count1, 5*d);
				}else if(gsl_vector_get(symOut,count1)>2*d&&gsl_vector_get(symOut, count1)<=4*d){
		gsl_vector_set(symOut, count1, 3*d);
				}else if(gsl_vector_get(symOut ,count1)>0&&gsl_vector_get(symOut, count1)<=2*d){
		gsl_vector_set(symOut, count1 , d);
				}

		}
	}
	//calculate mean square error
	gsl_blas_dcopy(symReceived, y_tmp);
	gsl_blas_dgemv(CblasNoTrans, 1 , pH, symOut, 0, symOut_tmp);
	gsl_vector_sub(y_tmp,symOut_tmp);
	*MSE=gsl_blas_dnrm2(y_tmp);


	gsl_matrix_free(G);
	gsl_permutation_free(p);
    gsl_vector_free(y_tmp);
    gsl_vector_free(symOut_tmp);
    gsl_matrix_free(LU);
    gsl_matrix_free(E);
	gsl_matrix_free(inverse);
    free(signum);
	return;
}



#endif /* MMSE_DETEC_H_ */
