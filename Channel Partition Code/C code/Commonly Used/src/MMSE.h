/*
 * MMSE.h
 *
 *  Created on: Dec 4, 2015
 *      Author: Preston Chen
 */
#include "RectangularQAMSlicer.h"
#ifndef MMSE_H_
#define MMSE_H_
void MMSE(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr, double pav, int M, gsl_vector_complex *psymOut)
{
gsl_complex alpha, beta1,beta2;
GSL_SET_COMPLEX(&alpha, 1,0);
GSL_SET_COMPLEX(&beta1, 1/snr, 0);
GSL_SET_COMPLEX(&beta2, 0, 0);
int Nr=pH->size1;
int Nt=pH->size2;
gsl_matrix_complex *G_pre, *G_preInv, *G;
gsl_permutation *p=gsl_permutation_calloc(Nt);
int signum[1];
signum[0]=1;
G_pre=gsl_matrix_complex_calloc(Nt, Nt);
G_preInv=gsl_matrix_complex_calloc(Nt, Nt);
G=gsl_matrix_complex_calloc(Nt,Nr);
gsl_matrix_complex_set_identity(G_pre);
gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH, pH, beta1, G_pre);
gsl_linalg_complex_LU_decomp(G_pre, p, signum);
gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);
gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, G_preInv, pH, beta2, G);
gsl_blas_zgemv(CblasNoTrans, alpha, G, preceived, beta2, psymOut);
RectangularQAMSlicer(psymOut, pav, M);
gsl_matrix_complex_free(G_pre);
gsl_matrix_complex_free(G_preInv);
gsl_matrix_complex_free(G);
gsl_permutation_free(p);



}




#endif /* MMSE_H_ */
