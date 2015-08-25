/*
 * SVR_DETECTOR.h
 *
 *  Created on: Jul 15, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SVR_DETECTOR_CPP_
#define SVR_DETECTOR_CPP_
#include "public.h"
#include"Initialization.h"
#include"WSS2_1Dsolver.h"
#include"TakeStep.h"
#include"Stopping_Criteria.h"
typedef struct{
	gsl_complex data;
	int index;

}Nalpha_Node;

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



void SVR_DETECTOR( gsl_vector_complex *preceived,  //received symbol vector
	 gsl_matrix_complex *pH, //channel propagation matrix
		float SNRb, // bit signal to noise ratio
		gsl_vector_complex *psymout, //detected transmitted symbol vector
		int start_M,  //initialization method
		int selection_M //work set selection strategy
		){
int count1, count2, count3; //counters

gsl_vector *alpha=gsl_vector_calloc(2*Nr);
gsl_vector *beta=gsl_vector_calloc(2*Nr);
gsl_matrix_complex *K_complex=gsl_matrix_complex_calloc(Nr,Nr);
gsl_blas_zherk(CblasUpper,CblasNoTrans, 1, pH,0, K_complex);  //calculate complex kernel matrix
gsl_matrix *K_r=gsl_matrix_calloc(Nr,Nr);   //real kernel matrix
gsl_matrix *K_i=gsl_matrix_calloc(Nr,Nr); //imaginary kernel matrix
/*
 * extract real and imaginary kernel matrix
 */
for(count1=0;count1<Nr;count1++){
	for(count2=0;count2<=count1;count2++){
		gsl_matrix_set(K_r,count1,count2,gsl_matrix_complex_get(K_complex,count1,count2).dat[0]);
		gsl_matrix_set(K_i,count1,count2,gsl_matrix_complex_get(K_complex,count1,count2).dat[1]);
	}
}

gsl_vector_complex *phi=gsl_vector_complex_calloc(Nr); //phi vector (intermediate parameter)
gsl_vector_complex *eta=gsl_vector_complex_calloc(Nr); //eta vector (intermediate parameter)
double L_pilot=0; //stopping criteria(L sub function)
double S_pilot=0; // stopping criteria(S sub function)
double G_pilot=0; //stopping criteria parameter

/*
 * Initialization
 */
//int M=Constellationsize;  //the size of modulation constellation
Initialization(alpha, beta,preceived, K_r, K_i, phi, eta, L_pilot, S_pilot, G_pilot, start_M);

/*
 * iterative update
 */
double L=0;
double H=C*(double)(1)/(double)(SNRb*log2(M));
int i_real, j_real, i_imag, j_imag;
count2=0; //iteration times of real part
count3=0;
int model_real=0;
int model_imag=1;
int label_real=0;
int label_imag=1;
double sigma1_gap_real;
double sigma1_gap_imag;
double sigma_sum_real;
double sigma2_gap_real;
double sigma2_gap_imag;
double sigma_sum_imag;
while(G_pilot>tol){

	i_real=j_real=0;
	sigma1_gap_real=0;
	sigma2_gap_imag=0;
	sigma_sum_real=0;
	WSS2_1Dsolver(alpha,phi,K_r,i_real,j_real,model_real);
	TakeStep(alpha,i_real,j_real, preceived, phi, eta, K_r, K_i, sigma1_gap_real, sigma2_gap_real, sigma_sum_real, label_real);


    i_imag=j_imag=0;
	sigma1_gap_imag=0;
	sigma2_gap_imag=0;
	sigma_sum_imag=0;
	WSS2_1Dsolver(beta, phi, K_r , i_imag, j_imag, model_imag);
	TakeStep(beta,i_imag ,j_imag ,preceived,phi,eta, K_r, K_i ,sigma1_gap_imag, sigma2_gap_imag,sigma_sum_imag,  label_imag);
    G_pilot=Stopping_Criteria(sigma1_gap_real, sigma1_gap_imag, sigma_sum_real, sigma2_gap_real, sigma2_gap_imag,
    		sigma_sum_imag, i_real, j_real, i_imag, j_imag,phi,eta, preceived, K_r, K_i, L_pilot, S_pilot);
    count3++; //record the iteration time

}

/*
 * reconstruction of detected symbol from dual variables
 */
gsl_vector_complex *psymout_temp1=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *psymout_temp2=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *alpha_complex=gsl_vector_complex_calloc(Nr);
gsl_vector_complex *beta_complex= gsl_vector_complex_calloc(Nr);
gsl_complex a_tmp;
gsl_complex b_tmp;
gsl_complex one;
gsl_complex zero;
gsl_complex one_i;
GSL_SET_COMPLEX(&one, 1, 0);
GSL_SET_COMPLEX(&zero, 0 ,0);
GSL_SET_COMPLEX(&one_i, 0, 1);
for(count2=0; count2<Nr;count2++){
	GSL_SET_COMPLEX(&a_tmp, gsl_vector_get(alpha,count2)-gsl_vector_get(alpha, count2+Nr), 0);
	gsl_vector_complex_set(alpha_complex, count2, a_tmp);
	GSL_SET_COMPLEX(&b_tmp, gsl_vector_get(beta,count2)-gsl_vector_get(beta, count2+Nr), 0);
    gsl_vector_complex_set(beta_complex, count2, b_tmp);
}
gsl_blas_zgemv(CblasConjTrans, one, pH, alpha_complex, zero, psymout_temp1);
gsl_blas_zgemv(CblasConjTrans, one, pH, beta_complex, zero, psymout_temp2);

gsl_blas_zaxpy(one_i, psymout_temp2, psymout_temp1);

/*
 * Rounding
 */
gsl_complex symout_temp;
double real_tmp;
double imag_tmp;
double d;  //constellation distance unit
for (count1=0;count1<Nt;count1++){
real_tmp=gsl_vector_complex_get(psymout_temp1,count1).dat[0];
imag_tmp=gsl_vector_complex_get(psymout_temp1,count1).dat[1];
d=sqrt(3/(2*Nt*(M-1)));
if(M==2){
	d=sqrt((double)(1)/(double)(Nt));
	if(real_tmp>0){
		GSL_SET_REAL(&symout_temp, d);
	}else{
		GSL_SET_REAL(&symout_temp, -d);
	}

   GSL_SET_IMAG(&symout_temp, 0);


}
else if(M==4){

	if(real_tmp>0){
		GSL_SET_REAL(&symout_temp, d);
	}else{
		GSL_SET_REAL(&symout_temp, -d);
	}

	if(imag_tmp>0){
		GSL_SET_IMAG(&symout_temp, d);
	}else{
		GSL_SET_IMAG(&symout_temp, -d);
	}

}
else if(M==16){
	if(real_tmp>0&&real_tmp<=2*d){
		GSL_SET_REAL(&symout_temp, d);
	}else if(real_tmp>2*d){
		GSL_SET_REAL(&symout_temp,3*d);
	}else if(real_tmp<=0&&real_tmp>-2*d){
		GSL_SET_REAL(&symout_temp,-d);
	}else if(real_tmp<=-2*d){
		GSL_SET_REAL(&symout_temp,-3*d);
	}

	if(imag_tmp>0&&imag_tmp<=2*d){
		GSL_SET_IMAG(&symout_temp, d);
	}else if(imag_tmp>2*d){
		GSL_SET_IMAG(&symout_temp,3*d);
	}else if(imag_tmp<=0&&real_tmp>-2*d){
		GSL_SET_IMAG(&symout_temp,-d);
	}else if(imag_tmp<=-2*d){
		GSL_SET_IMAG(&symout_temp,-3*d);
	}
}
else if(M==64){
	if(real_tmp>0&&real_tmp<=2*d){
		GSL_SET_REAL(&symout_temp, d);
	}else if(real_tmp>2*d&&real_tmp<=4*d){
		GSL_SET_REAL(&symout_temp,3*d);
	}else if(real_tmp>4*d&&real_tmp<=6*d){
		GSL_SET_REAL(&symout_temp,5*d);
	}else if(real_tmp>6*d){
		GSL_SET_REAL(&symout_temp,7*d);
	}else if(real_tmp<=0&&real_tmp>-2*d){
		GSL_SET_REAL(&symout_temp, -d);
	}else if(real_tmp<=-2*d&&real_tmp>-4*d){
		GSL_SET_REAL(&symout_temp,-3*d);
	}else if(real_tmp<=-4*d&&real_tmp>-6*d){
		GSL_SET_REAL(&symout_temp,-5*d);
	}else if(real_tmp<=-6*d){
		GSL_SET_REAL(&symout_temp,-7*d);
	}

	if(imag_tmp>0&&imag_tmp<=2*d){
		GSL_SET_IMAG(&symout_temp, d);
	}else if(imag_tmp>2*d&&imag_tmp<=4*d){
		GSL_SET_IMAG(&symout_temp,3*d);
	}else if(imag_tmp>4*d&&imag_tmp<=6*d){
		GSL_SET_IMAG(&symout_temp,5*d);
	}else if(imag_tmp>6*d){
		GSL_SET_IMAG(&symout_temp,7*d);
	}else if(imag_tmp<=0&&imag_tmp>-2*d){
		GSL_SET_IMAG(&symout_temp, -d);
	}else if(imag_tmp<=-2*d&&imag_tmp>-4*d){
		GSL_SET_IMAG(&symout_temp,-3*d);
	}else if(imag_tmp<=-4*d&&imag_tmp>-6*d){
		GSL_SET_IMAG(&symout_temp,-5*d);
	}else if(imag_tmp<=-6*d){
		GSL_SET_IMAG(&symout_temp,-7*d);
	}
}

gsl_vector_complex_set(psymout, count1, symout_temp);


}


/*
 * Free space
 */

gsl_vector_free(alpha);
gsl_vector_free(beta);
gsl_matrix_complex_free(K_complex);
gsl_matrix_free(K_r);
gsl_matrix_free(K_i);
gsl_vector_complex_free(phi);
gsl_vector_complex_free(eta);
gsl_vector_complex_free(psymout_temp1);
gsl_vector_complex_free(psymout_temp2);
gsl_vector_complex_free(alpha_complex);
gsl_vector_complex_free(beta_complex);

return;
}





#endif /* SVR_DETECTOR_CPP_ */
