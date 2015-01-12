unsigned long binaryerrors (unsigned long input, unsigned long output, unsigned long Q)
{
	unsigned long count=0UL, errors=0UL, temp=0UL, temp2=0UL;

	int M=0;
	M = (int) log2(static_cast<double>(Q));

	temp = input^output;
	errors = 0L;
	for (count=0; count<M; count++)
	  {
		temp2 = temp%2;
		temp = temp/2;
		if (temp2 != 0)
		 	errors = errors+1;

	  } 
	//printf ("\n input = %ld, \t output = %ld, \t errors = %ld", input, output, errors);

	return errors;
}
