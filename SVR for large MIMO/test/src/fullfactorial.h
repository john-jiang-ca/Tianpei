/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
/*
 * this fuction implement the full factorial design
 * 29/08/2014
 * Tianpei Chen
 */

#include "public.h"
void fullfact(
		int rho,  //the number of elements that use full expansion
//		int M,    //constellation size
		int *s_sub_index    //the indexes of the full expansion matrix  (pow(M,rho))
		)
{
	int ssize,ncycles,columns,nreps;
	int count1,count2,count3,count4;
	int *settings2,*settings1;
	settings2=(int*)calloc(1,int(pow(M,rho))*sizeof(int));
	ssize=int(pow(M,rho));
	ncycles=ssize;
	for(count1=0;count1<rho;count1++)
	{
//		memset(settings2,0,int(pow(M,rho))*sizeof(int));

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
//	memset(settings2,0,int(pow(M,rho))*sizeof(int));
	free(settings2);
//	free(settings1);



}
