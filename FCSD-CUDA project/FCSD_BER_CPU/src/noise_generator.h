void noise_generator (gsl_vector_complex *pnoise, gsl_rng *pr, double sigman)
{
	int Nr, count=0;
	Nr=pnoise->size;
	
	gsl_complex z;
	GSL_SET_COMPLEX (&z, 0.0, 0.0);

	for (count=0 ; count<Nr ; count++)
	    {
			GSL_SET_COMPLEX (&z, gsl_ran_gaussian (pr, sigman), gsl_ran_gaussian (pr, sigman));
			gsl_vector_complex_set (pnoise, count, z);
			//printf("z = %g \n", gsl_ran_gaussian (pr, 1.0));
	    }
	return;
		   
}
