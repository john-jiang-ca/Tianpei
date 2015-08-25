/*
 * Initialization.h
 *
 *  Created on: Aug 3, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 *      Initialization of SVR
 *      INPUT: alpha, alpha_hat, beta, beta_hat, preceived, method
 *      OUTPUT: phi, S_C_real, S_C_imag
 */
#include "common.h"
#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_
int Initialization(gsl_vector *alpha,//Lagrange multiplier
		gsl_vector *beta,//Lagrange multiplier
		const gsl_vector_complex *preceived,//received symbol vector(output of trainning data set)
		gsl_vector_complex *phi, // update parameter
		double S_C_real, //real stopping parameter
		double S_C_imag,  //imaginary stopping parameter
		int method //label which start strategy is used
		){
//	int method=0; //label of method
    int count1,count2;
    gsl_complex temp;
    GSL_SET_COMPLEX(&temp, 0, 0);
//cold start (L0-W0)
if(method==1){
for(count1=0;count1<Nr;count1++){
temp=gsl_vector_complex_get(preceived, count1);
gsl_vector_complex_set(phi, count1, temp);
S_C_real=0;
S_C_imag=0;  //need to consider what is stopping criteria of this SVR
}
}


	return method;
}




#endif /* INITIALIZATION_H_ */
