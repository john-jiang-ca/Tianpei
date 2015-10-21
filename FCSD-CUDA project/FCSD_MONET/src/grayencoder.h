void grayencoder (gsl_vector_ulong *pgraydata, gsl_vector_ulong *pgrayindexes, unsigned long Q)
{
	unsigned long count=0UL, count1=0UL, count2=0UL, temp=0UL;

	int M=0;
	M = (int) ceil(sqrt(Q));

	gsl_vector_ulong *ppamdata, *pQAMdata, *peqvdata;
	ppamdata = gsl_vector_ulong_calloc (M);
	pQAMdata = gsl_vector_ulong_calloc (Q);
	peqvdata = gsl_vector_ulong_calloc (Q);


	for (count=0; count<M; count++)
	  {
		count1=count;
		count2 = count;
		count2 = count2 >> 1;
		count1 = count1^count2;
		
		gsl_vector_ulong_set (ppamdata, count, count1);

	  }

	//vectorviewerint("ppamdata", ppamdata);
	
	for (count=0; count<Q; count++)
	  {
		count1 = M*gsl_vector_ulong_get(ppamdata, count/M) + gsl_vector_ulong_get (ppamdata, count%M);
		gsl_vector_ulong_set (pQAMdata, count, count1);
		gsl_vector_ulong_set (peqvdata, count1, count);   //the pdata is actually the gray code
		//according to the gray code we can find its indexs according to the indexes we can find the gray code

	  }

	//vectorviewerint ("pQAMdata", pQAMdata);
	//vectorviewerint ("peqvdata", peqvdata);
	gsl_vector_ulong_memcpy (pgraydata, peqvdata);
	gsl_vector_ulong_memcpy (pgrayindexes, pQAMdata);

	gsl_vector_ulong_free (ppamdata);
	gsl_vector_ulong_free (pQAMdata);
	gsl_vector_ulong_free (peqvdata);

	return;
}
