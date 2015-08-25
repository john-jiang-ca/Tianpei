#include"public.h"
void grayencoder (gsl_vector_ulong *pgraydata, gsl_vector_ulong *pgrayindexes, unsigned long Q)
{
	unsigned long count=0UL, value=0UL, count2=0UL, temp=0UL;

	int S=0;
	S = (int) ceil(sqrt(Q));

	gsl_vector_ulong *ppamdata, *pQAMdata, *peqvdata;
	ppamdata = gsl_vector_ulong_calloc (S);
	pQAMdata = gsl_vector_ulong_calloc (Q);
	peqvdata = gsl_vector_ulong_calloc (Q);
if(Q==4)
{
 for (count=0;count<Q;count++)
 {
	 switch (count) {
	 case 0:
	 gsl_vector_ulong_set(pQAMdata,count,0); break;
	 gsl_vector_ulong_set(peqvdata,0,count); break;
	 case 1:
		 gsl_vector_ulong_set(pQAMdata,count,1); break;
		 gsl_vector_ulong_set(peqvdata,1,count); break;
	 case 2:
		 gsl_vector_ulong_set(pQAMdata,count,3); break;
		 gsl_vector_ulong_set(peqvdata,3,count); break;
	 case 3:
		 gsl_vector_ulong_set(pQAMdata,count,2); break;
		 gsl_vector_ulong_set(peqvdata,2,count); break;
	 }
 }
}
else
{
	for (count=0; count<S; count++)
	  {
		value=count;
		count2 = count;
		count2 = count2 >> 1;
		value = value^count2;
		
		gsl_vector_ulong_set (ppamdata, count, value);

	  }

	//vectorviewerint("ppamdata", ppamdata);
	
	for (count=0; count<Q; count++)
	  {
		value = S*gsl_vector_ulong_get(ppamdata, count/S) + gsl_vector_ulong_get (ppamdata, count%S);
		gsl_vector_ulong_set (pQAMdata, count, value);
		gsl_vector_ulong_set (peqvdata, value, count);   //the pdata is actually the index
		//according to the gray code we can find its indexs according to the indexes we can find the gray code

	  }
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
