//============================================================================
// Name        : public.h
// Author      : Tianpei Chen
// Version     :
// Copyright   : Your copyright notice
// Description : Global static libraries and constant settings
//============================================================================

#ifndef PUBLIC_H_
#define PUBLIC_H_
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_ulong.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<time.h>
//#define MATRIX_SIZE 16
#define Nr 8
#define Nt 8
#define M  4//constellation size
#define IDC2D(i,j,row)    (i)*row+j
#define  C  1e-3     //Penalize weight for noise
#define epsilon 1e-5 //Training precision
#define tol 1e-4   //the tolerance for KKT condition
#define RUNCSVR_MIMO  //run CSVR-MIMO
#define RUNFCSD      //run FCSD
#define DEBUG  //debugging mode

#endif /* PUBLIC_H_ */
