/*
 * checkNonbound.h
 *
 *  Created on: Jul 16, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef CHECKNONBOUND_H_
#define CHECKNONBOUND_H_
#include "common.h"
Nalpha_Node checkNonbound(gsl_vector_complex *alpha){
//check the non bound alpha and store them in Nalpha
Nalpha_Node *Nalpha=(Nalpha_Node*)malloc(Constellationsize);
//  gsl_vector_complex *Nalpha=gsl_vector_complex_calloc(Constellationsize);
  return Nalpha;
}




#endif /* CHECKNONBOUND_H_ */
