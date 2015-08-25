#include "public.h"
void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, unsigned long Q)
{
//	int Nt_gener;
//	Nt_gener = pdata -> size;
	int count = 0;
	for (count=0; count<Nt ; count++)
	  {
		gsl_vector_ulong_set(pdata, count, gsl_rng_uniform_int (pr,Q));	
		/*printf ("\n data(%d)= %d \n", count, gsl_vector_ulong_get (pdata, count));*/	
	  }
	return;
}
