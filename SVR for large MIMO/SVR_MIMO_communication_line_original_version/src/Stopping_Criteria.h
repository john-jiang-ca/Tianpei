/*
 * Stopping_Criteria.h
 *
 *  Created on: Aug 20, 2015
 *      Author: Preston Chen
 *      This sub routine calculate duality gap
 */
#include "public.h"
#ifndef STOPPING_CRITERIA_H_
#define STOPPING_CRITERIA_H_
double Stopping_Criteria(
		double sigma1_gap_real, //sigma1-simga1_hat
		double sigma1_gap_imag, //tao1-tao1_hat
		double sigma_sum_real, //sigma1+sigma1_hat+sigma2+sigma2_hat
		double sigma2_gap_real, //sigma2-sigma2_hat
		double sigma2_gap_imag,//tao2-tao2_hat
		double sigma_sum_imag, //tao1+tao1_hat+tao2+tao2_hat
		int i_real, //index of sigma1
		int j_real, //index of sigma2
		int i_imag, //index of tao1
		int j_imag, //index of tao2
		gsl_vector_complex *phi, //intermediate vector
		gsl_vector_complex *eta, //intermediate vector
		gsl_vector_complex *preceived, //received symbol vector
		gsl_matrix *K_r, //real kernel matrix
		gsl_matrix *K_i, //imaginary kernel matrix
		double L, //L sub function
		double S //S sub function

		){
	int count1, count2, count3; //counters

 /*
 * update L sub funciton
 */
L+=sigma1_gap_real*(sigma1_gap_real*(gsl_matrix_get(K_r, i_real, i_real))-2*(gsl_vector_complex_get(phi, i_real).dat[0])
+gsl_vector_complex_get(preceived, i_real).dat[0]);

L+=sigma2_gap_real*(sigma2_gap_real*(gsl_matrix_get(K_r, j_real, j_real))-2*(gsl_vector_complex_get(phi, j_real).dat[0])
+gsl_vector_complex_get(preceived, j_real).dat[0])+2*sigma1_gap_real*sigma2_gap_real*gsl_matrix_get(K_r,i_real, j_real);

L+=sigma1_gap_imag*(sigma1_gap_imag*(gsl_matrix_get(K_r, i_imag, i_imag))-2*(gsl_vector_complex_get(phi, i_imag).dat[1])
+gsl_vector_complex_get(preceived, i_imag).dat[1]);

L+=sigma2_gap_imag*(sigma2_gap_imag*(gsl_matrix_get(K_r, j_imag, j_imag))-2*(gsl_vector_complex_get(phi, j_imag).dat[1])
+gsl_vector_complex_get(preceived, j_imag).dat[1])+2*sigma1_gap_imag*sigma2_gap_imag*gsl_matrix_get(K_r,i_imag, j_imag)+
		epsilon*(sigma_sum_imag+sigma_sum_real);
/*
 * update S sub function
 */
double S_temp;
for (count1=0; count1<Nr; count1++){
	S_temp=0;
	if(fabs(gsl_vector_complex_get(phi, count1).dat[0]+gsl_vector_complex_get(eta, count1).dat[0])>epsilon){
		S_temp=fabs(gsl_vector_complex_get(phi, count1).dat[0]+gsl_vector_complex_get(eta, count1).dat[0])-epsilon;
		S+=C*pow(S_temp, 2);
	}
	if(fabs(gsl_vector_complex_get(phi, count1).dat[1]+gsl_vector_complex_get(eta, count1).dat[1])>epsilon){
		S_temp=fabs(gsl_vector_complex_get(phi, count1).dat[1]+gsl_vector_complex_get(eta, count1).dat[1])-epsilon;
		S+=C*pow(S_temp, 2);
	}
}

 /*
  * update mutual term of duality gap
  */
/*what if we omit this term, that means we can calculate the approximate duality
 * without consider the imagineray part of kernel matrix
 */
S+=sigma1_gap_real*sigma1_gap_imag*gsl_matrix_get(K_i, i_real, i_imag)+sigma2_gap_real*sigma2_gap_imag
		*gsl_matrix_get(K_i, j_real, j_imag)+sigma1_gap_real*sigma2_gap_imag*gsl_matrix_get(K_i, i_real, j_imag)
+sigma2_gap_real*sigma1_gap_imag*gsl_matrix_get(K_i, j_real, i_imag)+sigma1_gap_real*gsl_vector_complex_get(eta, i_real).dat[1]
+sigma2_gap_real*gsl_vector_complex_get(eta, j_real).dat[1]+sigma1_gap_imag*(-gsl_vector_complex_get(eta, i_imag).dat[1])
+sigma2_gap_imag*(-gsl_vector_complex_get(eta, j_imag).dat[1]);


/*
 * update the global duality gap
 */

	return (L+S);
}




#endif /* STOPPING_CRITERIA_H_ */
