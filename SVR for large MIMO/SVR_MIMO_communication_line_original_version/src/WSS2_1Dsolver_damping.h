/*
 * WSS2_1Dsolver_damping.h
 *
 *  Created on: Aug 25, 2015
 *      Author: Preston Chen
 */
#include"public.h"

#ifndef WSS2_1DSOLVER_DAMPING_H_
#define WSS2_1DSOLVER_DAMPING_H_
void 	WSS2_1Dsolver_damping(gsl_vector *l_m,  //Lagrange Multiplier vector
        gsl_vector_complex *phi,  //intermediate update parameter
        gsl_matrix *kernel_r,    //real kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model //determine this routine is for real part (0) or imaginary part (1)
){

	double best_g, g_temp;
	double sigma, sigma_hat;
	double temp_l_m, temp_l_m_hat;
	gsl_vector_complex *phi_temp=gsl_vector_complex_calloc(Nr); //temporary intermediate phi vector
	double p_temp;
	gsl_blas_zcopy(phi, phi_temp);
	int count, count1;
	//double phi_p[2*Nr];
	//for(count=0;count<2*Nr;count++){  //take real and imaginary part of phi
	//	if(count<Nr){
	//  phi_p[count]=gsl_vector_complex_get(phi,count).dat[0];
	//	}else{
	//		phi_p[count]=gsl_vector_complex_get(phi,count-Nr).dat[1];
	//	}
	//}
	best_g=-1;
	F_i=S_i=0;
	sigma=0;
	sigma_hat=0;
	int first_round=0;
	int second_round=0;
	while(second_round==0){
	for(count=0;count<2*Nr;count++){
        if(first_round==1&&count==F_i){
     	   continue;
        }
		if(count<Nr){
			temp_l_m=gsl_vector_get(l_m, count)+(gsl_vector_complex_get(phi_temp,count).dat[model]-epsilon)/(gsl_matrix_get(kernel_r,count,count));
			if(temp_l_m<0){
				temp_l_m=0;
			}else if (temp_l_m>C){
				temp_l_m=C;
			}

			sigma=temp_l_m-gsl_vector_get(l_m,count);
			sigma_hat=0;
		}else{
			temp_l_m_hat=gsl_vector_get(l_m, count)-(double)((gsl_vector_complex_get(phi_temp,count).dat[model]+epsilon))/(double)((gsl_matrix_get(kernel_r, count,count)));
			if(temp_l_m_hat<0){
				temp_l_m_hat=0;
			}else if (temp_l_m_hat>C){
				temp_l_m_hat=C;
			}
			sigma=0;
			sigma_hat=temp_l_m_hat-gsl_vector_get(l_m,count);
		}

		g_temp=fabs((sigma-sigma_hat)*(-0.5*(sigma-sigma_hat)*gsl_matrix_get(kernel_r, count,count)+gsl_vector_complex_get(phi_temp, count).dat[model]-epsilon*(sigma+sigma_hat)/(sigma-sigma_hat)));
		if(g_temp>best_g){
			best_g=g_temp;
			if(first_round==0){
			F_i=count;
			}else{
				S_i=count;
			}
		}
	}

    if(first_round==0){
    	first_round=1;
    	/*
    	 * update phi intermediate vector
    	 */
for(count1=0;count1<Nr;count1++){
 p_temp=gsl_vector_complex_get(phi_temp, count1).dat[model]-(sigma-sigma_hat)*gsl_matrix_get(kernel_r,count1,F_i);
 gsl_vector_complex_get(phi_temp,count1).dat[model]=p_temp;
}

    }else{
    	second_round=1;
    	 p_temp=gsl_vector_complex_get(phi_temp, count1).dat[model]-(sigma-sigma_hat)*gsl_matrix_get(kernel_r,count1,S_i);
    	 gsl_vector_complex_get(phi_temp,count1).dat[model]=p_temp;

    }



	}

	gsl_blas_zcopy(phi_temp, phi);


}

#endif /* WSS2_1DSOLVER_DAMPING_H_ */
