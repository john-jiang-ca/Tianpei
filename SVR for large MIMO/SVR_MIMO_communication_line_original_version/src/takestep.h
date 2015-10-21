/*
 * takestep.h
 *
 *  Created on: Jul 16, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 *      this subroutine update 2 Lagrange variable pair sequentially
 */

#ifndef TAKESTEP_H_
#define TAKESTEP_H_
#include "public.h"
void TakeStep(gsl_vector *l_m,
		int index1, //coordinates that are chosen
		int index2, //coordinates that are chosen
		gsl_vector_complex *preceive, //receive symbol vector
		gsl_vector_complex *Phi, // phi complex vector
		gsl_vector_complex *eta, //eta complex vector
		gsl_matrix *K_r,// real kernel matrix
		gsl_matrix *K_i, //imaginary kernel matrix
//		double eta_sub, //mutual part
		double sigma1_gap,// the gap of the first coordinate pair
		double sigma2_gap,// the gap of the second coordinate pair
		double sigma_sum, //the sum of sigma
		int label, //0 real part 1 imaginary part
		double  L, //lower bound of dual variables
		double H //upper bound of dual variables
		){

	int count1, count2;  //counter
	int label_hat; //label determine whether to update alpha or alpha_hat
	double sigma1, sigma1_hat, sigma2, sigma2_hat;//step of dual variable
	double Phi_temp=0;
	double eta_temp=0;
	double l_m_temp=0;
	double S_sub_temp=0;
//update alpha or beta vector
/*
* update the first coordinate
*/
	if(index1<Nr){
		label_hat=0;
	}else{
		label_hat=1;
	}


l_m_temp=gsl_vector_get(l_m, index1)+pow(-1, label_hat)*((gsl_vector_complex_get(Phi, index1).dat[label]+pow(-1,label_hat+1)*epsilon)
		/gsl_matrix_get(K_r, index1, index1));

if(l_m_temp>C){
l_m_temp=C;
}else if(l_m_temp<0){
	l_m_temp=0;
}
if(label_hat==0 ){
	sigma1=gsl_vector_get(l_m,index1)-l_m_temp;
	sigma1_hat=0;
}else{
	sigma1=0;
	sigma1_hat=gsl_vector_get(l_m,index1)-l_m_temp;
}
sigma1_gap=sigma1-sigma1_hat;
gsl_vector_set(l_m, l_m_temp, index1);



/*
 * update the second coordinate
 */
if(index2<Nr){
	label_hat=0;
}else{
	label_hat=1;
}


l_m_temp=gsl_vector_get(l_m, index2)+pow(-1, label_hat)*((gsl_vector_complex_get(Phi, index2).dat[label]+pow(-1,label_hat+1)*epsilon)
	/gsl_matrix_get(K_r, index2, index2));

if(l_m_temp>C){
l_m_temp=C;
}else if(l_m_temp<0){
l_m_temp=0;
}
if(label_hat==0 ){
sigma2=gsl_vector_get(l_m,index2)-l_m_temp;
sigma2_hat=0;
}else{
sigma2=0;
sigma2_hat=gsl_vector_get(l_m,index2)-l_m_temp;
}
gsl_vector_set(l_m, l_m_temp, index2);
sigma2_gap=sigma2-sigma2_hat;


/*
 * update phi sub function
 * update eta sub function
 */
	for(count1=0;count1<Nr;count1++){
	Phi_temp-=(sigma1-sigma1_hat)*gsl_matrix_get(K_r, count1,index1 );
	Phi_temp-=(sigma2-sigma2_hat)*gsl_matrix_get(K_r, count1, index2);
	gsl_vector_complex_get(Phi, count1).dat[label]=Phi_temp;
    eta_temp+=pow(-1, label+1)*(sigma1-sigma1_hat)*gsl_matrix_get(K_r, count1, index1);
    eta_temp+=pow(-1, label+1)*(sigma2-sigma2_hat)*gsl_matrix_get(K_r, count1, index2);
    gsl_vector_complex_get(eta, count1).dat[label]=eta_temp;
	}


	sigma_sum=sigma1+sigma1_hat+sigma2+sigma2_hat;

	return ;
}


#endif /* TAKESTEP_H_ */
