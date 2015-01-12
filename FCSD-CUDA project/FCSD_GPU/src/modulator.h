
void modulator (gsl_vector_ulong *pdata, gsl_vector_complex *ptransmitted, gsl_vector_complex *psymbolconstellation)
{
	int Nt, count=0;
	Nt = pdata->size;
	
	for (count=0 ; count < Nt ; count++)
		gsl_vector_complex_set (ptransmitted, count, gsl_vector_complex_get (psymbolconstellation, gsl_vector_ulong_get (pdata, count)));
	

	return;
}

