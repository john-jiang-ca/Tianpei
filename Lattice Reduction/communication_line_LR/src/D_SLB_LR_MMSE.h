/*
 * Lattice Reduction (D_SLB) aided MMSE MIMO decoder
 * 4 QAM modulation test
 * INPUT
 * preceived: received symbol vector
 * pH: propagation matrix
 * SNR:signal to noise ratio in decimal
 * Nr: number of receive antennas
 * Nt: number of transmitt antennas
 * M: size of M-QAM modulation
 * psymbolconstellation: the symbol constellation
 * OUTPUT
 * symOut: the detected symbol vector 
 * Tianpei Chen
 * 2015-01-20
 */
#include"common.h"
#include<gsl/gsl_complex.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex_double.h>
//#include <gsl/gsl_matrix_complex.h>
//#include <gsl/gsl_vector_complex.h>
#define e 2.71828
void D_SLB_LR_MMSE(gsl_vector_complex *preceived,
		gsl_matrix_complex *pH,
		double SNR,
		int Nr,
		int Nt,
		int M,
		gsl_vector_complex *psymbolconstellation,
		gsl_vector_complex *symOut)
{
//Dual-Element Lattice Reduction-Smallest longest  Vector
int count1,count2,count3,count4;
int k;
gsl_complex one,zero,lamida,lamida_max;
//gsl_complex delta;
double delta, delta_temp,delta_temp1,delta_temp2,C_diag_max,C_diag_temp;
Nr=MATRIX_SIZE_Nr;
Nt=MATRIX_SIZE_Nt;
GSL_SET_COMPLEX(&one,1,0);
GSL_SET_COMPLEX(&zero,0,0);
GSL_SET_COMPLEX(&lamida,0,0);
//GSL_SET_COMPLEX(&delta,0,0);
gsl_matrix_complex *C=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *T=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *LU=gsl_matrix_complex_calloc(Nt,Nt);
gsl_vector_complex *C_diag=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *t_k=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *c_k=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *c_kk=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *t_i=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *c_i=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *c_ii=gsl_vector_complex_calloc(Nt);
gsl_matrix_complex_set_identity(T);
gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, one, pH, pH, one, LU);
gsl_permutation *p;
gsl_permutation_init(p);
int *x=(int*)calloc(1,sizeof(int));
*x=1;
gsl_linalg_complex_LU_decomp(LU,p,x);
gsl_linalg_complex_LU_invert(LU,p,C);
int *pilot=(int*)calloc(Nt,sizeof(int));
//pilot[k]=1;
bool determine=1;

do
{
for (count1=0;count1<Nt;count1++)
{
gsl_vector_complex_set(C_diag,count1,gsl_matrix_complex_get(C,count1,count1));
}
int i,k;
gsl_complex Ckk, Cik;
GSL_SET_COMPLEX(&Ckk,0,0);
GSL_SET_COMPLEX(&Cik,0,0);
//find the largest reducible Ckk
for(count2=0;count2<Nt;count2++)
{
	if(count2==1)
	{
		C_diag_max=0;
		k=1;
	}
	else
	{
		C_diag_temp=gsl_complex_abs(gsl_vector_complex_get(C_diag,count2));
		if(C_diag_temp>C_diag_max)
		{
			C_diag_max=C_diag_temp;
			k=count2;
		}
	}
}
Ckk=gsl_vector_complex_get(C_diag,k);
pilot[k]=0;
i=0;
for (count1=0;count1<Nt;count1++)
{
	if(count1==k)
	{
		continue;
	}
	
lamida.dat[0]=-round((gsl_complex_div(gsl_matrix_complex_get(C,count1,k),gsl_matrix_complex_get(C,count1,count1))).dat[0]);
lamida.dat[1]=-round((gsl_complex_div(gsl_matrix_complex_get(C,count1,k),gsl_matrix_complex_get(C,count1,count1))).dat[1]);
if(lamida.dat[0]==0&&lamida.dat[1]==0)
{
	pilot[count1]=1;
}
delta_temp1=gsl_complex_negative(gsl_complex_mul_real(gsl_matrix_complex_get(C,count1,count1),gsl_complex_abs2(lamida))).dat[0];
delta_temp2=gsl_complex_negative(gsl_complex_add(gsl_complex_mul(gsl_complex_conjugate(lamida),gsl_matrix_complex_get(C,count1,k)),gsl_complex_mul(lamida,gsl_matrix_complex_get(C,k,count1)))).dat[0];
delta=delta_temp1+delta_temp2;
if (count1==0)
{
	delta_temp=delta;
}
else
{
	if(delta>delta_temp)
	{
		delta_temp=delta;
		i=count1;
		lamida_max=lamida;
	}
}


}


//reduce the correlated columns
gsl_matrix_complex_get_col(t_k,T,k);
gsl_matrix_complex_get_col(t_i,T,i);
gsl_matrix_complex_get_col(c_k,C,k);
gsl_matrix_complex_get_col(c_i,C,i);
gsl_matrix_complex_get_row(c_kk,T,k);
gsl_matrix_complex_get_row(c_ii,T,i);
//gsl_vector_complex_scale(t_i,lamida);
//gsl_blas_zscal(lamida,t_i);
//gsl_blas_zscal(lamida,c_i);
//gsl_blas_zscale(gsl_complex_conjugate(lamida),c_ii);
gsl_blas_zaxpy(lamida,t_i,t_k);
gsl_blas_zaxpy(lamida,c_i,c_k);
gsl_blas_zaxpy(gsl_complex_conjugate(lamida),c_ii,c_kk);
gsl_matrix_complex_set_col(T,k,t_k);
gsl_matrix_complex_set_col(C,k,c_k);
gsl_matrix_complex_set_row(C,k,c_kk);
for(count1=0;count1<Nt;count1++)
{
 if(pilot[count1]==0)
 {
	 determine=0;
	 break;
 }
}

}
while(determine==0);
//update the matrix pH and T
gsl_matrix_complex_memcpy(LU,T);
gsl_linalg_complex_LU_decomp(LU,p,x);
gsl_linalg_complex_LU_invert(LU,p,T);
gsl_matrix_complex *I=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex_set_identity(I);
gsl_matrix_complex *T_update=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *pH_update=gsl_matrix_complex_calloc(Nr,Nt);
gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,one,T,I,zero,T_update );
gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,pH,T_update,zero,pH_update);

//release the storage
gsl_matrix_complex_free(C);
gsl_matrix_complex_free(T);
gsl_matrix_complex_free(LU);
gsl_vector_complex_free(C_diag);
gsl_vector_complex_free(t_k);
gsl_vector_complex_free(c_k);
gsl_vector_complex_free(c_kk);
gsl_vector_complex_free(t_i);
gsl_vector_complex_free(c_i);
gsl_vector_complex_free(c_ii);
gsl_matrix_complex_free(I);
gsl_matrix_complex_free(T_update);
gsl_matrix_complex_free(pH_update);


//MMSE linear detection 

//gsl_matrix_complex *I=gsl_matrix_calloc(Nt,Nt);
gsl_matrix_complex *wishart=gsl_matrix_complex_calloc(Nt,Nt);
//gsl_matrix_complex *LU=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *M_inverse=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex *G=gsl_matrix_complex_calloc(Nt,Nr);
gsl_vector_complex *x_hat=gsl_vector_complex_calloc(Nt);
gsl_vector_complex *x_hat_q=gsl_vector_complex_calloc(Nt);

gsl_matrix_complex_set_identity(I);
gsl_complex SNR_inv;
GSL_SET_COMPLEX(&SNR_inv,double(1/(SNR*log2(M))),0);
for  (count1=0;count1<Nt;count1++)
{
   gsl_matrix_complex_set(I,count1,count1,SNR_inv);
}

gsl_blas_zgemm(CblasConjTrans,CblasNoTrans, one, pH_update,pH_update,one,I);
gsl_linalg_complex_LU_decomp(I,p,x);
gsl_linalg_complex_LU_invert(I,p,M_inverse);
gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,one,M_inverse,pH_update,zero,G);
gsl_blas_zgemv(CblasNoTrans,one,G,preceived,zero,x_hat);
gsl_blas_zgemv(CblasNoTrans,one,T,x_hat,zero,x_hat_q);
//quantization
double temp;
double d=sqrt(3/(2*double(Nt)*(M-1)));
for (count1=0;count1<Nt;count1++)
{

  temp=gsl_vector_complex_get(x_hat_q,count1).dat[0];
  if(temp>0)
  {
          gsl_vector_complex_get(symOut,count1).dat[0]=d;
  }
  else
  {
	  gsl_vector_complex_get(symOut,count1).dat[0]=-d;
  }
  temp=gsl_vector_complex_get(x_hat_q,count1).dat[1];
  if(temp>0)
  {
	  gsl_vector_complex_get(symOut,count1).dat[1]=d;
  }
  else
  {
	  gsl_vector_complex_get(symOut,count1).dat[1]=-d;
  }


}

//free storage
gsl_matrix_complex_free(I);
gsl_matrix_complex_free(wishart);
gsl_matrix_complex_free(LU);
gsl_matrix_complex_free(M_inverse);
gsl_matrix_complex_free(G);
gsl_vector_complex_free(x_hat);
gsl_vector_complex_free(x_hat_q);
gsl_permutation_free(p);



}







