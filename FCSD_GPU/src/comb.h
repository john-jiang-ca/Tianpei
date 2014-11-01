/*
 * comb.h
 *
 *  Created on: Sep 10, 2014
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include"common.h"
#ifndef COMB_H_
#define COMB_H_
void comb(int m, int k, int row, int * aaa, int *subset, int *count3)
{
	  int count1,count2,count4;

//	  int *a;
	      for (count1=m;count1>=k;count1--)
	      {
	          aaa[k]=count1;
	          if (k>1)
	          {
	              comb(count1-1,k-1, row, aaa, subset, count3);
	          }
	          else
	          {
                   count4=0;
	              for (count2=aaa[0];count2>0;count2--)
	              {
	                  printf("%d ",aaa[count2]);
	                  subset[IDC2D(count3[0],count4,row)]=aaa[count2];
	                  count4++;


	               }
	              count3[0]++;
	              printf("\n");
	      }
	}
}



#endif /* COMB_H_ */
