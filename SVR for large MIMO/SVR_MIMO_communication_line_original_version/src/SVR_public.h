/*
 * SVR_public.h
 *
 *  Created on: Aug 21, 2015
 *      Author: Preston Chen
 */

#ifndef SVR_PUBLIC_H_
#define SVR_PUBLIC_H_
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
		);
void 	WSS2_1Dsolver(gsl_vector *l_m,  //Lagrange Multiplier vector
        gsl_vector_complex *phi,  //update parameter
        gsl_matrix *R_Kernel,    //real kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model //determine this routine is for real part (0) or imaginary part (1)
);

void 	WSS2_1Dsolver_damping(gsl_vector *l_m,  //Lagrange Multiplier vector
        gsl_vector_complex *phi,  //update parameter
        gsl_matrix *R_Kernel,    //real kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model //determine this routine is for real part (0) or imaginary part (1)
);
void TakeStep(gsl_vector *l_m,
		int index1, //coordinates that are chosen
		int index2, //coordinates that are chosen
		gsl_vector_complex *preceive, //receive symbol vector
		gsl_vector_complex *Phi, // phi complex vector
		gsl_vector_complex *eta, //eta complex vector
		gsl_matrix *K_r,// real kernel matrix
		gsl_matrix *K_i, //imaginary kernel matrix
		double sigma1_gap,// the gap of the first coordinate pair
		double sigma2_gap,// the gap of the second coordinate pair
		double sigma_sum, //the sum of sigma
		int label //0 real part 1 imaginary part
		);
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

		);




#endif /* SVR_PUBLIC_H_ */
