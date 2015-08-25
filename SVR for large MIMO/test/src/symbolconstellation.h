//============================================================================
// Name        : symbolconstellation.h
// Author      : Tianpei Chen
// Version     :
// Copyright   : Your copyright notice
// Description : Generate symbol constellation according to the index of the random generated data
// (pdata)
//============================================================================
#include "public.h"
void symbolconstellation(gsl_vector_complex *psymbolconstellation, int Nt_conste) {
	int Q;
	int count1=0;
	int count = 0;
	Q = psymbolconstellation->size;
	gsl_complex z;
	int S = 0;
	S = (int) ceil(sqrt(Q)); /*the square QAM of dimension Q is equivalent to 2 PAM of dimension S */
	/* The average energy of the square QAM is given by: Es=2/3*d^2*(Q-1)
	 where dmin=2d
	 Hence we get d=sqrt(3*Es/(2*(Q-1)))
	 We set Es=1/Nt_conste so the total transmitted power is normalized*/

	double d = 0.0;
	if (Q == 2)   //BPSK
			{
		d = sqrt((double)((double)(1) / (double)(Nt_conste)));
	}
	else {
		d = sqrt(3 / (2 * (float) (Nt_conste * (Q - 1))));
	}
	gsl_vector *psymbolconstellationreal;
	psymbolconstellationreal = gsl_vector_calloc(S);
		for (count = 0; count < S; count++) /*generates an S-PAM constellation*/
			gsl_vector_set(psymbolconstellationreal, count,
					d * (2 * count - S + 1));

		for (count = 0; count < Q; count++) {
			gsl_complex z;
			if (Q == 2) {
				GSL_SET_COMPLEX(&z,
						gsl_vector_get(psymbolconstellationreal, count), 0);
			} else {
				GSL_SET_COMPLEX(&z,
						gsl_vector_get(psymbolconstellationreal, count / S),
						gsl_vector_get(psymbolconstellationreal, (S-1)-count % S));
				gsl_vector_complex_set(psymbolconstellation, count, z);
			}
		}


	gsl_vector_free(psymbolconstellationreal);

	return;
}

