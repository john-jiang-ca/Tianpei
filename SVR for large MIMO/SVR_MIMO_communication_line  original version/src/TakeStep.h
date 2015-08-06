/*
 * TakeStep.h
 *
 *  Created on: Aug 5, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef TAKESTEP_H_
#define TAKESTEP_H_
#include"common.h"
void TakeStep(gsl_vector *l_m,
		gsl_vector_complex *phi,
		gsl_matrix *K,
		int label, // determine which part we are updating(real or imaginary)
		int i,
		int j
		){
	double l_m_temp, phi_temp;
	double sigmai=0,sigmai_hat=0, sigmaj=0, sigmaj_hat=0;
	gsl_complex phi_temp_complex;
	int count;
	if(i<Nr){
		l_m_temp=gsl_vector_get(l_m, i)+(gsl_vector_complex_get(phi,i).dat[label]-epsilon)
				/(gsl_matrix_get(K,i,i));
		if(l_m_temp<0){
			l_m_temp=0;
		}else if(l_m_temp>C){
			l_m_temp=C;
		}
		sigmai=gsl_vector_get(l_m,i)-l_m_temp;
		gsl_vector_set(l_m,i, l_m_temp);
	}else{
		l_m_temp=gsl_vector_get(l_m, i)-(gsl_vector_complex_get(phi,i).dat[label]+epsilon)
					/(gsl_matrix_get(K,i,i));
		if(l_m_temp<0){
			l_m_temp=0;
		}else if(l_m_temp>C){
			l_m_temp=C;
		}
		sigmai_hat=gsl_vector_get(l_m,i)-l_m_temp;
		gsl_vector_set(l_m,i, l_m_temp);
	}
	if(j<Nr){
		l_m_temp=gsl_vector_get(l_m, j)+(gsl_vector_complex_get(phi,j).dat[label]-epsilon)
				/(gsl_matrix_get(K,j,j));
		if(l_m_temp<0){
			l_m_temp=0;
		}else if(l_m_temp>C){
			l_m_temp=C;
		}
		sigmaj=gsl_vector_get(l_m,j)-l_m_temp;
		gsl_vector_set(l_m, j,l_m_temp);
	}else{
		l_m_temp=gsl_vector_get(l_m, j)-(gsl_vector_complex_get(phi,j).dat[label]-epsilon)
				/(gsl_matrix_get(K, j,j));
		if(l_m_temp<0){
			l_m_temp=0;
		}else if(l_m_temp>C){
			l_m_temp=C;
		}
		sigmaj_hat=gsl_vector_get(l_m,j)-l_m_temp;
		gsl_vector_set(l_m, j, l_m_temp);
	}
	/*
	 * updating phi vector stopping criteria
	 */
	for(count=0; count<Nr;count++){
      phi_temp=gsl_vector_complex_get(phi,count).dat[label]-(sigmai-sigmai_hat)*
     gsl_matrix_get(K, count, i)-(sigmaj-sigmaj_hat)*gsl_matrix_get(K,count, j);
      if(label==0){
    	  GSL_SET_COMPLEX(&phi_temp_complex, phi_temp,gsl_vector_complex_get(phi,count).dat[1]);

      }else{
    	 GSL_SET_COMPLEX(&phi_temp_complex, gsl_vector_complex_get(phi,count).dat[0],phi_temp);
      }
      gsl_vector_complex_set(phi, count ,phi_temp_complex);
	}
	return;
}




#endif /* TAKESTEP_H_ */
