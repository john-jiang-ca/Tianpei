/*
 * WSS2_1Dsolver.h
 *
 *  Created on: Aug 4, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 *      This is 2D selection subroutine based on maximum 2 1D direction
 *      subroutine, choose the two Lagrange multipliers with maximum 1D gain
 *      to perform 2D solver.
 *      INPUT: alpha(beta), alpha_hat(beta_hat), label1, label2, phi, kernel_r
 *      OUTPUT: First, Second
 */
#include "public.h"
#ifndef WSS2_1DSOLVER_H_
#define WSS2_1DSOLVER_H_
void 	WSS2_1Dsolver(gsl_vector *l_m,  //Lagrange Multiplier vector
        gsl_vector_complex *phi,  //intermediate vector
        gsl_vector_complex *eta, //intermediate vector
        gsl_matrix *kernel_r,    //real kernel matrix
        gsl_matrix *kernel_i, //imaginary kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model, //determine this routine is for real part (0) or imaginary part (1)
		gsl_vector_complex *phi_new, //new intermediate vector
		double sigma1, //dual variable1
		double sigma1_hat, //dual variable1_hat
		double sigma2, //dual variable2
		double sigma2_hat //dual variable2_hat
){
double best_g, g;
double sigma_tmp, sigma_hat_tmp;
double temp_l_m, temp_l_m_hat;
int count,count1;
double p_temp;
//gsl_vector_complex *phi_temp=gsl_vector_complex_calloc(Nr); //temporary intermediate phi vector
gsl_vector *dual_temp=gsl_vector_calloc(2*Nr); //temporary dual variable vector
gsl_vector *sig=gsl_vector_calloc(2*Nr);
//double phi_p[2*Nr];
//for(count=0;count<2*Nr;count++){  //take real and imaginary part of phi
//	if(count<Nr){
//  phi_p[count]=gsl_vector_complex_get(phi,count).dat[0];
//	}else{
//		phi_p[count]=gsl_vector_complex_get(phi,count-Nr).dat[1];
//	}
//}
best_g=-1;
F_i=S_i=-1;
sigma_tmp=0;
sigma_hat_tmp=0;
for(count=0;count<2*Nr;count++){
	if(count<Nr){
		temp_l_m=gsl_vector_get(l_m, count)+(gsl_vector_complex_get(phi,count).dat[model]-epsilon)/(gsl_matrix_get(kernel_r,count,count));
		if(temp_l_m<0){
			temp_l_m=0;
		}else if (temp_l_m>C){
			temp_l_m=C;
		}

		sigma_tmp=temp_l_m-gsl_vector_get(l_m,count);
		sigma_hat_tmp=0;
		 gsl_vector_set(dual_temp, count, temp_l_m);
        gsl_vector_set(sig, count, sigma_tmp);
	}else{
		temp_l_m_hat=gsl_vector_get(l_m, count)-(double)((gsl_vector_complex_get(phi,count-Nr).dat[model]+epsilon))/(double)((gsl_matrix_get(kernel_r, count-Nr,count-Nr)));
		if(temp_l_m_hat<0){
			temp_l_m_hat=0;
		}else if (temp_l_m_hat>C){
			temp_l_m_hat=C;
		}
		sigma_tmp=0;
		sigma_hat_tmp=temp_l_m_hat-gsl_vector_get(l_m,count);
		 gsl_vector_set(dual_temp, count, temp_l_m);
        gsl_vector_set(sig, count, sigma_hat_tmp);
	}
	int index_tmp;
   if(count<Nr){
	   index_tmp=count;
   }else{
	   index_tmp=count-Nr;
   }
	g=fabs((sigma_tmp-sigma_hat_tmp)*(-0.5*(sigma_tmp-sigma_hat_tmp)*gsl_matrix_get(kernel_r, index_tmp,index_tmp)
			+gsl_vector_complex_get(phi, index_tmp).dat[model]-epsilon*(sigma_tmp+sigma_hat_tmp)/(sigma_tmp-sigma_hat_tmp)));

	if(g>best_g){
		best_g=g;
		if(count!=F_i+Nr){
		S_i=F_i;
		}
		F_i=count;

	}
}

//update phi_new
double sigma_diff1=0;
double sigma_diff2=0;
int F_tmp=0;
int S_tmp=0;
if(F_i<Nr){
	sigma_diff1=gsl_vector_get(sig, F_i);
	F_tmp=F_i;
}else{
	sigma_diff1=-gsl_vector_get(sig, F_i);
	F_tmp=F_i-Nr;
}
if(S_i<Nr){
	sigma_diff2=gsl_vector_get(sig, S_i);
	S_tmp=S_i;
}else{
	sigma_diff2=-gsl_vector_get(sig, S_i);
	S_tmp=S_i-Nr;
}
for(count1=0;count1<Nr;count1++){

	p_temp=gsl_vector_complex_get(phi, count1).dat[model]-sigma_diff1*gsl_matrix_get(kernel_r,count1, F_tmp)
	-sigma_diff2*gsl_matrix_get(kernel_r, count1, S_tmp);
	gsl_vector_complex_get(phi_new, count1).dat[model]=p_temp;
}
gsl_vector_set(l_m, F_i, gsl_vector_get(dual_temp,F_i));
gsl_vector_set(l_m, S_i, gsl_vector_get(dual_temp,S_i));
if(F_i>Nr){
sigma1_hat=gsl_vector_get(sig,F_i);
sigma1=0;
F_i-=Nr;
}else{
	sigma1=gsl_vector_get(sig,F_i);
	sigma1_hat=0;
}

if(S_i>Nr){
sigma2_hat=gsl_vector_get(sig,S_i);
sigma2=0;
S_i-=Nr;
}else{
	sigma2=gsl_vector_get(sig,S_i);
	sigma2_hat=0;
}
/*
 * update eta
 */
double eta_temp;
for(count1=0;count1<Nr;count1++){
eta_temp+=pow(-1,model+1)*((sigma1-sigma1_hat)*gsl_matrix_get(kernel_i,count1,F_i)
	+(sigma2-sigma2_hat)*gsl_matrix_get(kernel_i,count1,S_i));
gsl_vector_complex_get(eta, count1).dat[(model^1)]=eta_temp;
}
//gsl_blas_zcopy(phi_new, phi_new);
//gsl_vector_complex_free(phi_new);
gsl_vector_free(dual_temp);
gsl_vector_free(sig);


}




#endif /* WSS2_1DSOLVER_H_ */
