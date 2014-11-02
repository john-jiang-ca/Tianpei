/*
 * fullfact.cuh
 *
 *  Created on: Oct 1, 2014
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "common.h"
#include<stdlib.h>
#include <stdio.h>
#include<string.h>
#include<math.h>
#include"common.h"
__device__ void fullfact_GPU(
		int rho,  //the number of elements that use full expansion
				int M,    //constellation size
				int *s_sub_index    //the indexes of the full expansion matrix  (pow(M,rho))
				)
{
	int ssize,ncycles,columns,nreps;
		int count1,count2,count3,count4;
		int *settings2;
		settings2=(int*)calloc(1,int(pow(M,rho))*sizeof(int));
		ssize=int(pow(M,rho));
		ncycles=ssize;
		for(count1=0;count1<rho;count1++)
		{
			memset(settings2,0,int(pow(M,rho))*sizeof(int));
		   int *settings1;
			nreps=int(ssize/ncycles);
			ncycles=ncycles/M;
			settings1=(int*)calloc(1,M*nreps*sizeof(int));
			for(count2=0;count2<M;count2++)
			{
				for(count3=0;count3<nreps;count3++)
				{
					settings1[count2*nreps+count3]=count2;
				}
			}


		for(count2=0;count2<ncycles;count2++)
		{
			for(count3=0;count3<M*nreps;count3++)
			{
				settings2[count2*M*nreps+count3]=settings1[count3];
			}
		}
		for(count4=0;count4<pow(M,rho);count4++)
		{
		s_sub_index[IDC2D(count1,count4,int(pow(M,rho)))]=settings2[count4];
		}

	    free(settings1);
		}
		memset(settings2,0,int(pow(M,rho))*sizeof(int));
		free(settings2);
}



#endif /* FULLFACT_CUH_ */
