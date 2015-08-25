#include "public.h"
void modulator (gsl_vector_ulong *pdata, gsl_vector_complex *ptransmitted, gsl_vector_complex *psymbolconstellation)
{
	int Nt_modu;
	int count=0;
	Nt_modu = pdata->size;
	
	for (count=0 ; count < Nt_modu ; count++)
		gsl_vector_complex_set (ptransmitted, count, gsl_vector_complex_get (psymbolconstellation, gsl_vector_ulong_get (pdata, count)));
	

	return;
}

