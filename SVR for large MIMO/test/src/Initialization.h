/*
 * Initialization.h
 *
 *  Created on: Aug 3, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 *      Initialization of SVR
 *      INPUT: alpha, alpha_hat, beta, beta_hat, preceived, method
 *      OUTPUT: phi, S_C_real, S_C_imag
 */
#include "public.h"
#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_
void  Initialization(gsl_vector *alpha,//Lagrange multiplier
		gsl_vector *beta,//Lagrange multiplier
		const gsl_vector_complex *preceived,//received symbol vector(output of trainning data set)
		gsl_matrix *K_r, //real kernel matrix
		gsl_matrix *K_i, //imaginary kernel matrix
		gsl_vector_complex *phi, // intermediate parameter
		gsl_vector_complex *eta, //intermediate parameter
		double L, //L sub function
		double S,  //S sub function
		double G_pilot, //global stopping parameter
		int method //label which start strategy is used
		){
//	int method=0; //label of method
    int count1,count2;
    gsl_complex temp;
    GSL_SET_COMPLEX(&temp, 0, 0);
//cold start (L0-W0)
if(method==1){
gsl_complex zero;
GSL_SET_COMPLEX(&zero, 0 ,0);

for(count1=0;count1<Nr;count1++){
temp=gsl_vector_complex_get(preceived, count1);
gsl_vector_complex_set(phi, count1, temp);
gsl_vector_complex_set(eta, count1, zero);
}
L=0;
S=0;
G_pilot=0;
}


	return;
}




#endif /* INITIALIZATION_H_ */
