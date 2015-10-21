/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
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
#include"fullfactorial.h"
#include"common.h"
void FCSD_CPU(
		gsl_vector_complex *preceived,   //received signal vector
		gsl_matrix_complex *pH,        //propagation matrix
				int Nt,             //number of transmit antennas
				int Nr,              //number of receive antennas
				int M,               //modulation scheme
		gsl_vector_complex *psymbolconstellation, //the symbol constellation
				float SNR,          //signal to noise ratio in decimal
		gsl_vector_complex *symOut,    //output symbol vector
		float *durationKernel
		)
{
	int count1, count2,count3, count4;
	int *list;
	list=(int*)calloc(1,MATRIX_SIZE*sizeof(int));
	for (count1=0;count1<Nt;count1++)
	{
		list[count1]=count1+1;
	}
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
//		    printf("the pH is:\n");
//		    for (count1=0;count1<MATRIX_SIZE;count1++)
//			    {
//			    	for(count2=0;count2<MATRIX_SIZE;count2++)
//			    	{
//			    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pH,count1,count2).dat[0],gsl_matrix_complex_get(pH,count1,count2).dat[1]);
//			    	}
//			    	printf("\n");
//			    }
//
//		    printf("\n");
        	gsl_matrix_complex *pH_permuted,*pH_loop;    //permuted propagation channel matrix
		    pH_permuted=gsl_matrix_complex_calloc(MATRIX_SIZE,MATRIX_SIZE);
		    pH_loop=gsl_matrix_complex_calloc(MATRIX_SIZE,MATRIX_SIZE);
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
//						    		pH_loop[IDC2D(count1,count3,MATRIX_SIZE)]=pH[IDC2D(count1,count3,MATRIX_SIZE)];
						    		gsl_matrix_complex_set(pH_loop,count1,count3,gsl_matrix_complex_get(pH,count1,count3));
						    	}
						    }

					}
//					printf("the pH_loop is:\n");
//				    for (count1=0;count1<MATRIX_SIZE;count1++)
//					    {
//					    	for(count3=0;count3<MATRIX_SIZE-count2;count3++)
//					    	{
////						    		pH_loop[IDC2D(count1,count3,MATRIX_SIZE)]=pH[IDC2D(count1,count3,MATRIX_SIZE)];
//					    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pH_loop,count1,count3).dat[0],gsl_matrix_complex_get(pH_loop,count1,count3).dat[1]);
//					    	}
//					    	printf("\n");
//					    }
                gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,alpha,pH_loop,pH_loop,beta,Pprod );
//			    printf("the wishart matrix is:\n");
		//        MATRIX_COLUMNROWT(Pprod,Pprod_column,MATRIX_SIZE-count2);
//			    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
//			    {
//			    	for(mm=0;mm<MATRIX_SIZE-count2;mm++)
//			    	{
//			    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(Pprod,count1,mm).dat[0],gsl_matrix_complex_get(Pprod,count1,mm).dat[1]);
//			    	}
//			    	printf("\n");
//			    }
//			    printf("the first element of wishart matrix is:\n");
//			    printf("%0.4f%+0.4fi\n",gsl_matrix_complex_get(Pprod,0,0).dat[0],gsl_matrix_complex_get(Pprod,0,0).dat[1]);
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
//			      printf("the duration of matrix inverse CPU is %0.4f:\n", duration);
//			      printf("%0.4f ", CLOCKS_PER_SEC);
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
//		    printf("the inverse matrix of wishart matrix is\n");
//			    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
//			    {
//			    	for(mm=0;mm<MATRIX_SIZE-count2;mm++)
//			    	{
//			    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(Pinv,count1,mm).dat[0],gsl_matrix_complex_get(Pinv,count1,mm).dat[1]);
//			    	}
//			    	printf("\n");
//			    }
//			    printf("the first element of inverse wishart matrix is:\n");
//			    	    printf("%0.4f%+0.4fi\n",gsl_matrix_complex_get(Pinv,0,0).dat[0],gsl_matrix_complex_get(Pinv,0,0).dat[1]);
		//	    MATRIX_ROWCOLUMNT(Pinv,Pinv_temp);

			    for (count1=0;count1<MATRIX_SIZE-count2;count1++)
			    {

			    		gsl_vector_complex_set(diag,count1,gsl_matrix_complex_get(Pinv,count1,count1));

			    }
//			    printf("the diagonal elements of the inverse wishart matrix is\n");
			    int m;
//			    for (m=0;m<MATRIX_SIZE-count2;m++)
//			    {
//			    	printf("%0.4f+i%0.4f ", gsl_vector_complex_get(diag,m).dat[0],gsl_vector_complex_get(diag,m).dat[1]);
//			    }
//		        printf("this is for test\n");
               for(count1=0;count1<MATRIX_SIZE-count2;count1++)
               {
            	   temp=sqrt(pow(gsl_vector_complex_get(diag,count1).dat[0],2)+pow(gsl_vector_complex_get(diag,count1).dat[1],2));
            	   gsl_vector_set(diag_abs,count1,temp);
               }
//			    printf("the diagonal elements are:\n");
//			    for(mm=0;mm<MATRIX_SIZE-count2;mm++)
//			    {
//			    	printf("%f ", gsl_vector_get(diag_abs,mm));
//			    }
//			    printf("\n");
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
//			    printf("the new list is:\n");
//			    for(mm=0;mm<MATRIX_SIZE;mm++)
//			    {
//			    	printf("%d ", list[mm]);
//			    }
//			    printf("\n");
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

//		       printf("\n");
//		       printf("the renewed version of the propagation matrix is:\n");
//			    for (count1=0;count1<MATRIX_SIZE;count1++)
//			    {
//			    	for(mm=0;mm<MATRIX_SIZE-count2-1;mm++)
//			    	{
//			    		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pH_loop,count1,mm).dat[0],gsl_matrix_complex_get(pH_loop,count1,mm).dat[1]);
//			    	}
//			    	printf("\n");
//			    }

		  gsl_matrix_complex_free(Pprod);
		  gsl_matrix_complex_free(Pinv);
		  gsl_matrix_complex_free(pH_temp1);
          gsl_vector_complex_free(diag);
          gsl_vector_free(diag_abs);
          gsl_matrix_complex_free(LU);

			}
			 gsl_matrix_complex_free(pH_loop);
			end=clock();
			duration=double(start-end);
//update the permuted matrix
//printf("the ordering list of CPU is:\n");
//for(count1=0;count1<Nt;count1++)
//{
//	printf("%d ", list[count1]);
//}
//printf("\n");

			 free(list_temp);
		for(count1=0;count1<MATRIX_SIZE;count1++)
		{
		for(count2=0;count2<MATRIX_SIZE;count2++)
		{
		gsl_matrix_complex_set(pH_permuted,count2,(MATRIX_SIZE-count1-1),gsl_matrix_complex_get(pH,count2,(list[count1]-1)));
		}
		}
//QR or cholesky factorization
gsl_matrix_complex *pW=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *pR=gsl_matrix_complex_calloc(Nr,Nt);
gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,alpha,pH_permuted,pH_permuted,beta,pW);
gsl_matrix_complex_memcpy(pR,pW);
double durationchol;
//clock_t start,end;
start=clock();
gsl_linalg_complex_cholesky_decomp(pR);
end=clock();
durationchol=double((end-start));
//printf("the duration of chol CPU is %0.4f ms:\n", durationchol);
for(count1=0;count1<Nt;count1++)
{
	for(count2=0;count2<count1;count2++)
	{
		gsl_matrix_complex_set(pR,count1,count2,beta);
	}
}
//printf("the upper triangular matrix is:\n");
//for(count1=0;count1<Nt;count1++)
//{
//	for(count2=0;count2<Nt;count2++)
//	{
//		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(pR,count1,count2).dat[0],gsl_matrix_complex_get(pR,count1,count2).dat[1]);
//	}
//	printf("\n");
//}



//get the unconstrained estimation s_hat
gsl_vector_complex *s_hat=gsl_vector_complex_calloc(Nt);
gsl_matrix_complex *I=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *I_temp=gsl_matrix_complex_calloc(Nt,Nt);
//gsl_matrix_complex *LU1=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *I_inv=gsl_matrix_complex_calloc(Nt,Nt);
gsl_permutation *f=gsl_permutation_calloc(Nt);
gsl_complex snr;
GSL_SET_COMPLEX(&snr,1/SNR,0);

	for(count1=0;count1<Nt;count1++)
	{
		gsl_matrix_complex_set(I,count1,count1,snr);
	}
//	printf("the matrix I get 1 is:\n");
//	for(count1=0;count1<Nt;count1++)
//	{
//	for( count2=0;count2<Nt;count2++)
//	{
//		printf("%0.4f%+0.4fi ", gsl_matrix_complex_get(I,count2,count1).dat[0],gsl_matrix_complex_get(I,count2,count1).dat[1]);
//	}
//	printf("\n");
//	}
gsl_matrix_complex_add(I,pW);
gsl_permutation_init(f);
gsl_linalg_complex_LU_decomp(I,f,x);
gsl_linalg_complex_LU_invert(I,f, I_inv);
gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,alpha,I_inv,pH_permuted,beta,I_temp);
gsl_blas_zgemv(CblasNoTrans,alpha,I_temp,preceived,beta,s_hat);
gsl_matrix_complex_free(pW);
gsl_permutation_free(f);
gsl_matrix_complex_free(I_temp);
gsl_matrix_complex_free(I);
//gsl_matrix_complex_free(LU1);
gsl_matrix_complex_free(I_inv);
//printf("the s_hat of CPU is:\n");
int count=0;
//for( count=0;count<Nt;count++)
//{
//	printf("%0.4f%+0.4fi ", gsl_vector_complex_get(s_hat,count).dat[0],gsl_vector_complex_get(s_hat,count).dat[1]);
//}
//printf("\n");
gsl_matrix_complex_free(pH_permuted);

//ful expansion and single expansion
		int rho=ceil(sqrt(Nt)-1);
		int pathNum=pow(M,rho);
		float d;    //the minimum distance of the signal constellation is 2d
//		gsl_vector *Eu=gsl_vector_calloc(pathNum);
		int Eu;
		float Eu_mini;
		float Eu_value;
//		gsl_matrix_complex *s_matrix_Q=gsl_matrix_complex_calloc(Nt,pathNum);
		int *s_sub_index=(int*)calloc(1,rho*pow(M,rho)*sizeof(int));
//		memset(s_sub_index,0,rho*pow(M,rho)*sizeof(int));
			fullfact(rho,M,s_sub_index);    //get  the indexes of all the possible rho length symbol vectors
start=clock();
gsl_vector *Eu_vector=gsl_vector_calloc(pathNum);
for(count1=0;count1<pathNum;count1++)
{
	gsl_vector_complex *s_sub=gsl_vector_complex_calloc(rho);
	gsl_vector_complex *s_Q=gsl_vector_complex_calloc(Nt);
	gsl_vector_complex *Eu_temp=gsl_vector_complex_calloc(Nt);
	gsl_vector_complex *s_matrix_Q=gsl_vector_complex_calloc(Nt);
	gsl_complex S_temp;
	for(count3=0;count3<rho;count3++)
	{
//	s_sub_share[IDC2D((tx),count1,Nt)]=psymbolconstellation[s_sub_index[IDC2D((tx),count1,Nt)]];
		gsl_vector_complex_set(s_sub,count3,gsl_vector_complex_get(psymbolconstellation,s_sub_index[IDC2D(count3,count1,pathNum)]));
	}
//	printf("the full expansion sub symbol vector is:\n");
//	for(count4=0;count4<rho;count4++)
//	{
//		printf("%0.4f%+0.4fi ", gsl_vector_complex_get(s_sub,count4).dat[0],gsl_vector_complex_get(s_sub,count4).dat[1]);
//	}
//	printf("\n");
	for(count2=(Nt-1);count2>=0;count2--)
	{
		if (count2<Nt-rho)
				{
					GSL_SET_COMPLEX(&S_temp,gsl_vector_complex_get(s_hat,count2).dat[0],gsl_vector_complex_get(s_hat,count2).dat[1]);
					for (count3=count2+1;count3<Nt; count3++)
					{
						S_temp=gsl_complex_add(S_temp,gsl_complex_mul(gsl_complex_div(gsl_matrix_complex_get(pR,count2,count3),gsl_matrix_complex_get(pR,count2,count2)),(gsl_complex_sub(gsl_vector_complex_get(s_hat,count3),gsl_vector_complex_get(s_matrix_Q,count3)))));
					}
//					s_share[IDC2D(threadIdx.x,count1,Nt)]=R_temp[threadIdx.x];
					// Quantization according to corresponding signal constellation
							if(M==2)   //BPSK
								{

							    	 d=sqrt(double(double(1)/double(Nt)));
							    	 gsl_complex BPSK_temp;
							    	 if(S_temp.dat[0]>0)
							    	 {
							    		  BPSK_temp.dat[0]=d;
							    		  BPSK_temp.dat[1]=0;
							    	 }
							    	 else
							    	 {
							    		  BPSK_temp.dat[0]=-d;
							    		  BPSK_temp.dat[1]=0;
							    	 }
					//		    	 GSL_SET_COMPLEX(&BPSK_temp,float(int(S_temp.dat[0]/d)*d)+float(int(round(S_temp.dat[0]/d-int(S_temp.dat[0]/d)))*d),0);
					//		    	 gsl_matrixsqrt(3/(2* (float)(Nt*(M-1))));_complex_get(s_matrix_Q,count2,pathNum).dat[0]=float(int(S_temp.dat[0]/d)*d)+float(int(round(S_temp.dat[0]/d-int(S_temp.dat[0]/d)))*d);
					//		    	 gsl_matrix_complex_get(s_matrix_Q,count2,pathNum).dat[1]=0;
							    	 gsl_vector_complex_set(s_matrix_Q,count2,BPSK_temp);

								}
								else if(M==4)   //QPSK
								{
									d=sqrt(3/(2* (float)(Nt*(M-1))));
					               gsl_complex QAM4;
									if(S_temp.dat[0]<0)
					                {
					                 GSL_SET_REAL(&QAM4,-d);
					                }
									else
									{
										GSL_SET_REAL(&QAM4,d);
									}
									if(S_temp.dat[1]<0)
					                {
					                 GSL_SET_IMAG(&QAM4,-d);
					                }
									else
									{
										GSL_SET_IMAG(&QAM4,d);
									}
									gsl_vector_complex_set(s_matrix_Q,count2,QAM4);

								}

								else if(M==16)  //16QAM
								{
									gsl_complex QAM16_temp;
								  d=sqrt(3/(2* (float)(Nt*(M-1))));
					//			  GSL_SET_COMPLEX(&QAM16_temp,2*d*float(int(S_temp.dat[0]/(2*d)))+d*(S_temp.dat[1]/sqrt(pow(S_temp.dat[0],2))),2*d*int(S_temp.dat[1]/(2*d))+d*(S_temp.dat[1]/sqrt(pow(S_temp.dat[1],2))));
					//			  GSL_SET_COMPLEX(&QAM16_temp,float(int((S_temp.dat[0]-d)/(2*d))*2*d+d)+float(round((S_temp.dat[0]-d)/(2*d)-int((S_temp.dat[0]-d)/(2*d))))*2*d,float(int((S_temp.dat[1]-d)/(2*d))*2*d+d)+float(round((S_temp.dat[1]-d)/(2*d)-int((S_temp.dat[1]-d)/(2*d))))*2*d);
					//			  if(S_temp.dat[0]<=(-3*d))
					//			  {
					//				QAM16_temp.dat[0]=-3*d;
					//			  }
					//			  else if(S_temp.dat[0]>=3*d)
					//			  {
					//				  QAM16_temp.dat[0]=3*d;
					//			  }
					//
					//			  if(S_temp.dat[1]<=(-3*d))
					//			  {
					//				QAM16_temp.dat[1]=-3*d;
					//			  }
					//			  else if(S_temp.dat[1]>=3*d)
					//			  {
					//				  QAM16_temp.dat[1]=3*d;
					//			  }



								  			  if(S_temp.dat[0]<(-2*d))
								  			  {
								  				QAM16_temp.dat[0]=(-3*d);
								  			  }
								  			  else if(S_temp.dat[0]>(2*d))
								  			  {
								  				  QAM16_temp.dat[0]=(3*d);
								  			  }
								  			  else if(S_temp.dat[0]>=0&&S_temp.dat[0]<=2*d)
								  			  {
								  				QAM16_temp.dat[0]=d;
								  			  }
								  			  else if(S_temp.dat[0]>=(-2*d)&&S_temp.dat[0]<=0)
								  			  {
								  				QAM16_temp.dat[0]=(-d);
								  			  }




								  			  if(S_temp.dat[1]<(-2*d))
								  			  {
								  				QAM16_temp.dat[1]=(-3*d);
								  			  }
								  			  else if(S_temp.dat[1]>(2*d))
								  			  {
								  				  QAM16_temp.dat[1]=(3*d);
								  			  }
								  			  else if(S_temp.dat[1]>=0&&S_temp.dat[1]<=(2*d))
								  			  {
								  				QAM16_temp.dat[1]=d;
								  			  }
								  			  else if(S_temp.dat[1]>=(-2*d)&&S_temp.dat[1]<=0)
								  			  {
								  				QAM16_temp.dat[1]=(-d);
								  			  }
					//			  			  printf("my solution 1 is:\n");
					//			  			  printf("%0.4f%+0.4fi ", QAM16_temp1.dat[0],QAM16_temp1.dat[1]);
					//			  			  printf("my solution 2 is:\n");
					//			  			  printf("%0.4f%+0.4fi ", QAM16_temp.dat[0],QAM16_temp.dat[1]);
					              gsl_vector_complex_set(s_matrix_Q,count2,QAM16_temp);
								}
								else if(M==64)   //64QAM
									{
					                   gsl_complex QAM64_temp;
								    	 d=sqrt(3/(2* (float)(Nt*(M-1))));
								    	 //real part
								    	 if(S_temp.dat[0]<=(-6*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,-7*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(-4*d)&&S_temp.dat[0]>(-6*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,-5*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(-2*d)&&S_temp.dat[0]>(-4*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,-3*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(0)&&S_temp.dat[0]>(-2*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,-d);
								    	 }
								    	 if(S_temp.dat[0]>(6*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,7*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(6*d)&&S_temp.dat[0]>(4*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,5*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(4*d)&&S_temp.dat[0]>(2*d))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,3*d);
								    	 }
								    	 else if(S_temp.dat[0]<=(2*d)&&S_temp.dat[0]>(0))
								    	 {
								    		 GSL_SET_REAL(&QAM64_temp,d);
								    	 }
								    	 //image part
								    	 if(S_temp.dat[1]<=(-6*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,-7*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(-4*d)&&S_temp.dat[1]>(-6*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,-5*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(-2*d)&&S_temp.dat[1]>(-4*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,-3*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(0)&&S_temp.dat[1]>(-2*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,-d);
								    	 }
								    	 if(S_temp.dat[1]>(6*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,7*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(6*d)&&S_temp.dat[1]>(4*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,5*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(4*d)&&S_temp.dat[1]>(2*d))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,3*d);
								    	 }
								    	 else if(S_temp.dat[1]<=(2*d)&&S_temp.dat[1]>(0))
								    	 {
								    		 GSL_SET_IMAG(&QAM64_temp,d);
								    	 }
//								    	 GSL_SET_COMPLEX(&QAM64_temp,float(int(S_temp.dat[0]/d)*d)+float(int(round(S_temp.dat[0]/d-int(S_temp.dat[0]/d)))*d),float(int(S_temp.dat[1]/d)*d)+float(int(round(S_temp.dat[1]/d-int(S_temp.dat[1]/d)))*d));
								    	 gsl_vector_complex_set(s_matrix_Q,count2,QAM64_temp);
								     }
				}
				else
				{
//					S_temp=s_sub[(Nt-count1-1)];
					GSL_SET_COMPLEX(&S_temp,gsl_vector_complex_get(s_sub,(Nt-count2-1)).dat[0],gsl_vector_complex_get(s_sub,(Nt-count2-1)).dat[1]);
					gsl_vector_complex_set(s_matrix_Q,count2,S_temp);
				}
}
	//calculating the Euclidean distance
	for(count3=0;count3<Nt;count3++)
		{
			gsl_vector_complex_set(s_Q,count3,gsl_vector_complex_get(s_matrix_Q,count3));
		}
		gsl_vector_complex_sub(s_Q,s_hat);
		gsl_blas_zgemv(CblasNoTrans,alpha,pR,s_Q,beta,Eu_temp);
//		gsl_vector_set(Eu,count1,gsl_blas_dznrm2(Eu_temp));
		Eu_value=gsl_blas_dznrm2(Eu_temp);

		gsl_vector_set(Eu_vector,count1,Eu_value);
if(count1==0)
{
	Eu_mini=Eu_value;
	Eu=0;
	 for(count3=0;count3<Nt;count3++)
	 	{
	 		gsl_vector_complex_set(symOut,(list[count3]-1),gsl_vector_complex_get(s_matrix_Q,(Nt-count3-1)));
	 	}
}
else
{
 if(Eu_value<Eu_mini)
 {
	 Eu_mini=Eu_value;
	 Eu=count1;
	 for(count3=0;count3<Nt;count3++)
	 	{
	 		gsl_vector_complex_set(symOut,(list[count3]-1),gsl_vector_complex_get(s_matrix_Q,(Nt-count3-1)));
	 	}

 }
}




        gsl_vector_complex_free(s_matrix_Q);
		gsl_vector_complex_free(s_sub);
		gsl_vector_complex_free(s_Q);
		gsl_vector_complex_free(Eu_temp);
		s_matrix_Q=NULL;
        s_sub=NULL;
        s_Q=NULL;
        Eu_temp=NULL;
}
int mini;
mini=gsl_vector_min_index(Eu_vector);
gsl_vector_free(Eu_vector);
//printf("the original Eu is %d:\n", Eu);
//printf("the minimum Eu is %d:\n", mini);
//printf("the Eu value is:\n");
//for(count1=0;count1<pathNum;count1++)
//{
//	printf("%f ", gsl_vector_get(Eu_vector,count1));
//}
//printf("\n");
end=clock();
*durationKernel=float(end-start);







free(list);
list=NULL;
free(x);
x=NULL;
gsl_matrix_complex_free(pR);
pR=NULL;
gsl_vector_complex_free(s_hat);
s_hat=NULL;
//gsl_matrix_complex_free(I);
//gsl_matrix_complex_free(s_matrix_Q);
//gsl_vector_free(Eu);
free(s_sub_index);
s_sub_index=NULL;

}




