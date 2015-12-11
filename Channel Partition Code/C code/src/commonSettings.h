/*
 * commonSettings.h
 *
 *  Created on: Dec 4, 2015
 *      Author: Preston Chen
 */

#ifndef COMMONSETTINGS_H_
#define COMMONSETTINGS_H_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix_complex_float.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_ulong.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "data_generator.h"
#include "grayencoder.h"
#include "symbolconstellation.h"
#include "modulator.h"
#include "channel_generator.h"
#include "corr_matrix_generator.h"
#include "error_channel_generator.h"
#include "noise_generator.h"
#include "MMSE.h"
#include "symErrorCheck.h"
#include "demodulator.h"
#include "binaryerrors.h"
#include "RectangularQAMSlicer.h"
#include "MMSE_OSIC.h"
//#define DEBUG    //debugging mode
#define fileName "/home/tchen44/code/CommunicationTest/test data/MMSE_OSIC_test.txt"  //the output file
int Corr_Ind=0; //the correlation channel mode (0 close 1 open)
int Est_Ind=0; //channel estimation error mode (0 close 1 open)
int receiveAntennas=32;     //number of receive antenna
int transmitAntennas=32;     //number of transmit antenna
int symConstellationSize=16;     //modulation scheme
double Rr=0;   // receive correlation parameter. Should be a real number in the interval between 0 and 1
double Rt=0;	// transmit correlation parameter. Should be a real number in the interval between 0 and 1
double gammasq=0; //the parameter for the channel estimation error
int minSymErrors=200; //the minimum symbol error accumulated
int minChannelRealizations=1e4; //the minimum channel realizations
int Start_SNR=0;   //the start SNR   (receive SNR in dB)
int End_SNR=16;     //the end SNR
int Step_SNR=2;   //the step of SNR


#endif /* COMMONSETTINGS_H_ */
