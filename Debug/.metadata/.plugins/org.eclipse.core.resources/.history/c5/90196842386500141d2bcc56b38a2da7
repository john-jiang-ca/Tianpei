/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include"common.h"
#include"cu_complex_operation.cuh"
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include<cublas_v2.h>
#include<cuda_runtime.h>
#include<cuda.h>
#include<cuda_profiler_api.h>
#include<cudaProfiler.h>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_ulong.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_permutation.h>
void FCSD_ordering_CPU(
		cuComplex *pH,
		int *list,
		cuComplex *d_pH_permuted
		)
{
	int count1, count2,count3, count4;
//		int *list;
//		list=(int*)calloc(1,MATRIX_SIZE*sizeof(int));

		//FCSD ordering
	    int *x=(int*)calloc(1,sizeof(int));
	    *x=1;

	//    gsl_matrix_complex *pH_temp;
	//			pH_temp=gsl_matrix_complex_calloc(Nr,Nt);
				int j,i,mm;
				int *list_temp;
				list_temp=(int*)calloc(1,MATRIX_SIZE*sizeof(int));
	//			j=(int*)malloc(sizeof(int));
				int N=ceil(sqrt(MATRIX_SIZE)-1);
	//			for( count=1;count<MATRIX_SIZE+1;count++)
	//			{
	//				list[count-1]=count;
	//			}
//			    printf("the pH is:\n");
//			    for (count1=0;count1<MATRIX_SIZE;count1++)
//				    {
//				    	for(count2=0;count2<MATRIX_SIZE;count2++)
//				    	{
//				    		printf("%0.4f%+0.4fi ", pH[IDC2D(count1,count2,MATRIX_SIZE)].x,pH[IDC2D(count1,count2,MATRIX_SIZE)].y);
//				    	}
//				    	printf("\n");
//				    }
//
//			    printf("\n");
			    gsl_matrix_complex *pH_loop;
	            pH_loop=gsl_matrix_complex_calloc(MATRIX_SIZE,MATRIX_SIZE);
	        	gsl_matrix_complex *pH_permuted;    //permuted propagation channel matrix
			    pH_permuted=gsl_matrix_complex_calloc(MATRIX_SIZE,MATRIX_SIZE);
			    gsl_complex alpha,beta;
			    float temp;
			    GSL_SET_COMPLEX(&alpha,1,0);
	            GSL_SET_COMPLEX(&beta,0,0);
	            clock_t start,end;
	            double duration;
	            start=clock();
				for(count2=0;count2<MATRIX_SIZE-1;count2++)
				{
					gsl_matrix_complex *Pprod,*Pinv,*pH_temp1,*LU;
					gsl_vector_complex *diag;
					gsl_vector *diag_abs;
					gsl_complex a;
					Pprod=gsl_matrix_complex_calloc((MATRIX_SIZE-count2),(MATRIX_SIZE-count2));
					Pinv=gsl_matrix_complex_calloc((MATRIX_SIZE-count2),(MATRIX_SIZE-count2));
					pH_temp1=gsl_matrix_complex_calloc(MATRIX_SIZE,(MATRIX_SIZE-count2-1));
					diag=gsl_vector_complex_calloc((MATRIX_SIZE-count2));
					diag_abs=gsl_vector_calloc((MATRIX_SIZE-count2));
				    LU=gsl_matrix_complex_calloc((MATRIX_SIZE-count2),(MATRIX_SIZE-count2));
	//		    memset(Pprod,0, MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));

						if(count2==0)
						{
						    for (count1=0;count1<MATRIX_SIZE;count1++)
							    {
							    	for(count3=0;count3<MATRIX_SIZE;count3++)
							    	{
							    		GSL_SET_COMPLEX(&a,pH[IDC2D(count1,count3,MATRIX_SIZE)].x,pH[IDC2D(count1,count3,MATRIX_SIZE)].y);
	//						    		pH_loop[IDC2D(count1,count3,MATRIX_SIZE)]=pH[IDC2D(count1,count3,MATRIX_SIZE)];
							    		gsl_matrix_complex_set(pH_loop,count1,count3,a);
							    	}
							    }

						}
//						printf("the pH_loop is:\n");
//					    for (count1=0;count1<MATRIX_SIZE;count1++)
//						    {
//						    	for(count3=0;count3<MATRIX_SIZE-count2;count3++)
//						    	{
//	//						    		pH_loop[IDC2D(count1,count3,MATRIX_SIZE)]=pH[IDC2D(count1,count3,MATRIX_SIZE)];
//						    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pH_loop,count1,count3).dat[0],gsl_matrix_complex_get(pH_loop,count1,count3).dat[1]);
//						    	}
//						    	printf("\n");
//						    }
	                gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,alpha,pH_loop,pH_loop,beta,Pprod );
//				    printf("the wishart matrix is:\n");
			//        MATRIX_COLUMNROWT(Pprod,Pprod_column,MATRIX_SIZE-count2);
//				    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
//				    {
//				    	for(mm=0;mm<MATRIX_SIZE-count2;mm++)
//				    	{
//				    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(Pprod,count1,mm).dat[0],gsl_matrix_complex_get(Pprod,count1,mm).dat[1]);
//				    	}
//				    	printf("\n");
//				    }
//				    printf("the first element of wishart matrix is:\n");
//				    printf("%0.4f%+0.4fi\n",gsl_matrix_complex_get(Pprod,0,0).dat[0],gsl_matrix_complex_get(Pprod,0,0).dat[1]);
				    //to be completed by matrix inverse function
				    int column=MATRIX_SIZE-count2;
				    int row=MATRIX_SIZE;
			//	    (cuComplex*)realloc(Pinv_temp,(MATRIX_SIZE-count2)*(MATRIX_SIZE-count2)*sizeof(cuComplex));
			//	    (cuComplex*)realloc(Pprod_temp1,(MATRIX_SIZE-count2)*(MATRIX_SIZE-count2)*sizeof(cuComplex));
			//	    MATRIX_COLUMNROWT(Pprod,Pprod_temp1);
			//	    Pprod=Pprod_temp1;
	//			    printf("the Pprod is\n");
	//			    for (count1=0;count1<MATRIX_SIZE;count1++)
	//			 	    {
	//			 	    	for(mm=0;mm<MATRIX_SIZE;mm++)
	//			 	    	{
	//			 	    		printf("%0.4f%+0.4fi ", Pprod[IDC2D(count1,mm,MATRIX_SIZE)].x,Pprod[IDC2D(count1,mm,MATRIX_SIZE)].y);
	//			 	    	}
	//			 	    	printf("\n");
	//			 	    }
	//			    MATRIX_INVERSE(Pprod,Pinv,row,column);
				    gsl_permutation *p=gsl_permutation_calloc((MATRIX_SIZE-count2));

				    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
					    {
					    	for(count3=0;count3<MATRIX_SIZE-count2;count3++)
					    	{
	//						    		pH_loop[IDC2D(count1,count3,MATRIX_SIZE)]=pH[IDC2D(count1,count3,MATRIX_SIZE)];
					    		gsl_matrix_complex_set(LU,count1,count3,gsl_matrix_complex_get(Pprod,count1,count3));
					    	}
					    }
				    gsl_permutation_init(p);
				    double duration;
				    clock_t start, end;
				    start=clock();
				    gsl_linalg_complex_LU_decomp(LU,p,x);
				    gsl_linalg_complex_LU_invert(LU,p, Pinv);
				    end=clock();
				    duration=double((end-start));
//				      printf("the duration of matrix inverse CPU is %0.4f:\n", duration);
//				      printf("%0.4f ", CLOCKS_PER_SEC);
	//			    gsl_matrix_complex *I;
	//			    I=gsl_matrix_complex_calloc(MATRIX_SIZE,MATRIX_SIZE);
	//			    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,alpha,Pprod,Pinv,beta, I);
	//			    printf("the Identity matrix is\n");
	//					    for (count1=0;count1<MATRIX_SIZE;count1++)
	//					    {
	//					    	for(mm=0;mm<MATRIX_SIZE;mm++)
	//					    	{
	//					    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(I,count1,mm).dat[0],gsl_matrix_complex_get(I,count1,mm).dat[1]);
	//					    	}
	//					    	printf("\n");
	//					    }
				    gsl_permutation_free(p);
//			    printf("the inverse matrix of wishart matrix is\n");
//				    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
//				    {
//				    	for(mm=0;mm<MATRIX_SIZE-count2;mm++)
//				    	{
//				    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(Pinv,count1,mm).dat[0],gsl_matrix_complex_get(Pinv,count1,mm).dat[1]);
//				    	}
//				    	printf("\n");
//				    }
//				    printf("the first element of inverse wishart matrix is:\n");
//				    	    printf("%0.4f%+0.4fi\n",gsl_matrix_complex_get(Pinv,0,0).dat[0],gsl_matrix_complex_get(Pinv,0,0).dat[1]);
			//	    MATRIX_ROWCOLUMNT(Pinv,Pinv_temp);

				    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
				    {

				    		gsl_vector_complex_set(diag,count1,gsl_matrix_complex_get(Pinv,count1,count1));

				    }
//				    printf("the diagonal elements of the inverse wishart matrix is\n");
//				    int m;
//				    for (m=0;m<MATRIX_SIZE-count2;m++)
//				    {
//				    	printf("%0.4f+i%0.4f ", gsl_vector_complex_get(diag,m).dat[0],gsl_vector_complex_get(diag,m).dat[1]);
//				    }
//			        printf("this is for test\n");
	               for(count1=0;count1<MATRIX_SIZE-count2;count1++)
	               {
	            	   temp=sqrt(pow(gsl_vector_complex_get(diag,count1).dat[0],2)+pow(gsl_vector_complex_get(diag,count1).dat[1],2));
	            	   gsl_vector_set(diag_abs,count1,temp);
	               }

					if(count2<N)
					{
						j=gsl_vector_max_index(diag_abs);    //full expansion
					}
					else
					{
						j=gsl_vector_min_index(diag_abs);              //single expansion
					}

				    //record the index of the symbol vector
				    list_temp[count2]=list[count2+j];   //choose the *j th element in the new list to be detected first
				    list[count2+j]=0;
				for (count3=0;count3<MATRIX_SIZE-count2-1;count3++)
				{
					for(int count4=count2;count4<MATRIX_SIZE;count4++)
					{
						if(list[count4]!=0)
						{
					list_temp[count2+count3+1]=list[count4];   //choose the next MATRIX_SIZE-(count2+1) to form the next sequence, the sequence order is the same
					list[count4]=0;
						break;
						}
					}
				}
				for(count3=0;count3<MATRIX_SIZE;count3++)
				{
				list[count3]=list_temp[count3];
				}
//				    printf("the new list is:\n");
//				    for(mm=0;mm<MATRIX_SIZE;mm++)
//				    {
//				    	printf("%d ", list[mm]);
//				    }
			//update the propagation matrix by remove the jth column in pH_temp

			//	 pH_temp1_temp=
			//	    pH_temp1=(cuComplex*)realloc(pH_temp1,(MATRIX_SIZE-count2-1)*MATRIX_SIZE*sizeof(cuComplex));
			//	 if(pH_temp1_temp!=NULL)
			//	 {
			//	 pH_temp1=pH_temp1_temp;
			//	 }
			//		else
			//	    	{
			//	    		printf("realloc for pH_temp1_temp failed\n");
			//	    	}
				for(count3=0;count3<MATRIX_SIZE;count3++)
				 {
					 for(count1=0;count1<MATRIX_SIZE-count2;count1++)
					 {
						 if(count1<(j))
						 {
	//		              pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]=pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
							 gsl_matrix_complex_set(pH_temp1,count3,count1,gsl_matrix_complex_get(pH_loop,count3,count1));
						 }
						 else if(count1>(j))
						 {
	//						 pH_temp1[IDC2D(count3,count1-1,(MATRIX_SIZE-count2-1))]=pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2))];
							 gsl_matrix_complex_set(pH_temp1,count3,(count1-1),gsl_matrix_complex_get(pH_loop,count3,count1));
						 }
					 }
				 }
				gsl_matrix_complex_free(pH_loop);
				pH_loop=gsl_matrix_complex_calloc(MATRIX_SIZE,(MATRIX_SIZE-count2-1));
				for(count3=0;count3<MATRIX_SIZE;count3++)
				 {
					 for(count1=0;count1<MATRIX_SIZE-count2-1;count1++)
					 {
	//		         pH_loop[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))]= pH_temp1[IDC2D(count3,count1,(MATRIX_SIZE-count2-1))];
						 gsl_matrix_complex_set(pH_loop,count3,count1,gsl_matrix_complex_get(pH_temp1,count3,count1));
					 }
				 }

//			       printf("\n");
//			       printf("the renewed version of the propagation matrix is:\n");
//				    for (count1=0;count1<MATRIX_SIZE;count1++)
//				    {
//				    	for(mm=0;mm<MATRIX_SIZE-count2-1;mm++)
//				    	{
//				    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pH_loop,count1,mm).dat[0],gsl_matrix_complex_get(pH_loop,count1,mm).dat[1]);
//				    	}
//				    	printf("\n");
//				    }
			  gsl_matrix_complex_free(Pprod);
			  gsl_matrix_complex_free(Pinv);
			  gsl_matrix_complex_free(pH_temp1);
	          gsl_vector_complex_free(diag);
	          gsl_vector_free(diag_abs);
	          gsl_matrix_complex_free(LU);

				}
				end=clock();
				duration=double(start-end);
	//update the permuted matrix
	cuComplex *pH_permuted1;
	pH_permuted1=(cuComplex*)malloc(MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex));
	for(count1=0;count1<MATRIX_SIZE;count1++)
	{
	for(count2=0;count2<MATRIX_SIZE;count2++)
	{
	pH_permuted1[IDC2D(count1,(MATRIX_SIZE-count2-1),MATRIX_SIZE)]=pH[IDC2D(count1,list[count2]-1,MATRIX_SIZE)];
	}
	}
  cudaMemcpy(d_pH_permuted,pH_permuted1,MATRIX_SIZE*MATRIX_SIZE*sizeof(cuComplex),cudaMemcpyHostToDevice);
  free(pH_permuted1);
  gsl_matrix_complex_free(pH_loop);
  gsl_matrix_complex_free(pH_permuted);
  free(list_temp);
  free(x);
//    ret=cublasDestroy(handle);
//    if (ret != CUBLAS_STATUS_SUCCESS)
//       {
//           printf("cublasDestory returned error code %d, line(%d)\n", ret, __LINE__);
//           exit(EXIT_FAILURE);
//       }
}


