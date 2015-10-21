void constrcandidates (gsl_complex sbabai, double bound, double d, int M, gsl_matrix *prings, gsl_matrix_int *pringsindexes, gsl_vector_int *pringsize, gsl_vector *pringradius, int numberofrings, gsl_vector_complex *psymbolconstellation, int maxcandidates, double alpha, gsl_complex rll, int *pfoundcandidates, gsl_vector_complex *pcandidates, gsl_vector_ulong *poutputcandidates, gsl_vector *pdistances, int *ptotalcandidates)
{
	int Q, count=0, ringsize=0, count2=0, numberofcandidates=0, temp=0, ii=0;
	Q = psymbolconstellation -> size;

	double Rc=0.0, Rest=0.0, Rcsq=0.0, Restsq=0.0, thetac=0.0, thetaest=0.0, Rllsq=0.0, alphasq=0.0, weight=0.0, bias=0.0, temp2=0.0, etha=0.0, Rmax=0.0, Rmaxsq=0.0;
	Rest = gsl_complex_abs (sbabai);
	Restsq = pow(Rest, 2.0);

	thetaest = gsl_complex_arg (sbabai);

	Rllsq = gsl_complex_abs2(rll);

	alphasq = pow (alpha, 2.0);
	Rmax = gsl_vector_get (pringradius, numberofrings-1);
	Rmaxsq = pow (Rmax, 2.0);
	//printf("Rmaxsq = %g \n", Rmaxsq);

	//printf ("bound = %g \n", bound);
	//printf ("sbabai = %g + (%g)i \n ", GSL_REAL(sbabai), GSL_IMAG(sbabai));

	gsl_vector_complex *pcurrentcandidates;
	pcurrentcandidates = gsl_vector_complex_calloc (Q);

	gsl_vector_ulong *pcurrentoutputs;
	pcurrentoutputs = gsl_vector_ulong_calloc (Q);
	
	gsl_vector *pclosest;
	pclosest = gsl_vector_calloc (Q);
	

	gsl_complex z, candidate;
		GSL_SET_COMPLEX (&z, 0.0, 0.0);
		GSL_SET_COMPLEX (&candidate, 0.0, 0.0);
	//z = gsl_vector_complex_get (psymbolconstellation, 0);
	//Rmax = gsl_complex_abs(z);
	//Rmaxsq = pow(Rmax, 2.0);
	//printf("Rmaxsq = %g \n", Rmaxsq);

	gsl_vector_view vtemp;

	long index1=0L, index2=0L;

	numberofcandidates = 0;
	//printf("numberofcandidates = %d \n", numberofcandidates);

	if (maxcandidates ==1)
	 {
		if (GSL_REAL (sbabai) <= d*(2-M))
			index1 = 0L;
		else
	  	{
			if (GSL_REAL (sbabai) >= d*(M-2))
				index1 = M-1L;
			else
				index1 = lrint ((double)((GSL_REAL (sbabai)-d*(1-M))/(2*d)));
	  	}
			   
		if (GSL_IMAG (sbabai) <= d*(2-M))
			index2 = 0L;
		else
	  	{
			if (GSL_IMAG (sbabai) >= d*(M-2))
				index2 = M-1L;
			else
				index2 = lrint ((double)((GSL_IMAG (sbabai)-d*(1-M))/(2*d)));
	  	}

		candidate = gsl_vector_complex_get (psymbolconstellation, (unsigned long)(index1*M+index2));
		Rc = gsl_complex_abs (candidate);
		Rcsq = pow(Rc, 2.0);
		bias = alphasq*(Rcsq-Rmaxsq);

		thetac = gsl_complex_arg (candidate);
		etha = 1.0/(2*Rc*Rest)*(Rcsq+Restsq-1.0/Rllsq*(bound+bias));

		if (cos(thetac-thetaest) >= etha)
		 {
			gsl_vector_complex_set (pcurrentcandidates, numberofcandidates, candidate);
			gsl_vector_ulong_set (pcurrentoutputs, numberofcandidates, (unsigned long) (index1*M+index2));
			z = gsl_complex_sub (sbabai, candidate);

			weight = Rllsq*gsl_complex_abs2(z)-bias;

			gsl_vector_set (pclosest, numberofcandidates, weight);
			numberofcandidates++;

		 }

	 }
	else
	 {
		for (count=0; count<numberofrings; count++)
	 	 {
			Rc = gsl_vector_get (pringradius, count);
			Rcsq = pow(Rc, 2.0);
			//printf("Rc = %g \n", Rc);

			bias = alphasq*(Rcsq-Rmaxsq);

			etha = 1.0/(2*Rc*Rest)*(Rcsq+Restsq-1.0/Rllsq*(bound+bias));

			//printf("etha = %g \n", etha);
			ringsize = gsl_vector_int_get (pringsize, count);

			if (etha > 1.0)
				continue;
			else 
		 	{
				if (etha < -1.0)
			 	 {
					for (count2=0; count2<ringsize; count2++)
				 	 {
						temp = gsl_matrix_int_get (pringsindexes, count, count2);
						z = gsl_vector_complex_get (psymbolconstellation, temp);
						gsl_vector_complex_set (pcurrentcandidates, numberofcandidates, z); 
						gsl_vector_ulong_set (pcurrentoutputs, numberofcandidates, (unsigned long) temp);
					
						z = gsl_complex_sub (sbabai, z);

						weight = Rllsq*gsl_complex_abs2(z)-bias;

						gsl_vector_set (pclosest, numberofcandidates, weight);
						numberofcandidates++;
				 	 }
				//numberofcandidates += count2+1;
			 	 }
				else 
			 	 {

					for (count2=0; count2 <ringsize; count2++)
				 	 {
						temp2 = gsl_matrix_get (prings, count, count2)-thetaest;
						//printf ("temp2 = %g \n", temp2);
						temp2 = cos(temp2);

						if (temp2 >= etha)
					 	 {
							temp = gsl_matrix_int_get (pringsindexes, count, count2);
							z = gsl_vector_complex_get (psymbolconstellation, temp);
							//printf("R = %g \n", gsl_complex_abs(z));
							gsl_vector_complex_set (pcurrentcandidates, numberofcandidates, z); 
							gsl_vector_ulong_set (pcurrentoutputs, numberofcandidates, (unsigned long) temp);
					
							z = gsl_complex_sub (sbabai, z);

							weight = Rllsq*gsl_complex_abs2(z)-bias;

							gsl_vector_set (pclosest, numberofcandidates, weight);
							numberofcandidates++;
					 	 }
				 	 }
					//numberofcandidates += ii;
				 }
			 }
		 }
	  }
	//printf("numberofcandidates = %d \n", numberofcandidates);
	//printf ("good = %d \n", numberofcandidates>0);
	//vectorviewer("pcurrentcandidates", pcurrentcandidates);
	//vectorviewerreal("pclosest", pclosest);

	if (numberofcandidates)
	 {
		gsl_permutation *p;
		vtemp = gsl_vector_subvector(pclosest, 0, numberofcandidates);
		p = gsl_permutation_calloc (numberofcandidates);
		gsl_sort_vector_index (p, &vtemp.vector);

		//vectorviewerint("p", p);

		for (ii=0; ii<numberofcandidates; ii++)
		 {
			count = gsl_permutation_get (p, ii);
			//printf("count=%d \n", count);
			gsl_vector_complex_set (pcandidates, ii, gsl_vector_complex_get(pcurrentcandidates, count));
			gsl_vector_ulong_set (poutputcandidates, ii, gsl_vector_ulong_get(pcurrentoutputs, count));
			gsl_vector_set (pdistances, ii, gsl_vector_get(pclosest, count));
		 }

		gsl_permutation_free (p);

		*pfoundcandidates = 1;
		*ptotalcandidates = numberofcandidates;
			
	 }
	else
	{
	  *pfoundcandidates = 0;
	  *ptotalcandidates = 0;
	}

	//printf ("good = %d \n", numberofcandidates>0);
	//vectorviewer("pcandidates", pcandidates);
	//vectorviewerreal("pdistances", pdistances);

	gsl_vector_complex_free (pcurrentcandidates);
	gsl_vector_ulong_free (pcurrentoutputs);
	gsl_vector_free (pclosest);



	return;
}
