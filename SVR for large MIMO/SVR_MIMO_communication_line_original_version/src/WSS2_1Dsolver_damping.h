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
        gsl_matrix *R_Kernel,    //real kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model //determine this routine is for real part (0) or imaginary part (1)
){


}

#endif /* WSS2_1DSOLVER_DAMPING_H_ */
