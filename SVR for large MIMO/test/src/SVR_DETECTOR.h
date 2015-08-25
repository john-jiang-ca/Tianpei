/*
 * SVR_DETECTOR.h
 *
 *  Created on: Aug 21, 2015
 *      Author: Preston Chen
 */
#include"public.h"
#ifndef SVR_DETECTOR_H_
#define SVR_DETECTOR_H_


void SVR_DETECTOR( gsl_vector_complex *preceived,  //received symbol vector
	 gsl_matrix_complex *pH, //channel propagation matrix
		float SNRb, // bit signal to noise ratio
		gsl_vector_complex *psymout, //detected transmitted symbol vector
		int start_M,  //initialization method
		int selection_M //work set selection strategy
		);


#endif /* SVR_DETECTOR_H_ */
