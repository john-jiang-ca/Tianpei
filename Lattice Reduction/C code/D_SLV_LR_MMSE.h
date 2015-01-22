/*
 * Lattice Reduction aided Belief propagation MIMO decoder
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
#define e 2.71828
D_SLB_LR_BP(preceived,pH,SNR,Nr,Nt,M,psymbolconstellation,symOut)
{
//Dual-Element Lattice Reduction-Smallest longest  Vector
int count1,count2,count3,count4;
gsl_complex one,zero,lamida;
//gsl_complex delta;
double delta, delta_temp;
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
gsl_vector_complex *t_i=gsl_vector_complex_calloc(Nt)
gsl_vector_complex *c_i=gsl_vector_complex_calloc(Nt)
gsl_vector_complex *c_ii=gsl_vector_complex_calloc(Nt)
gsl_matrix_complex_set_identity(T);
gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, one, pH, pH, one, LU);
gsl_permutation *p;
gsl_permutation_init(p);
int *x=(int*)calloc(1,sizeof(int));
*x=1;
gsl_linalg_complex_LU_decomp(LU,p,x);
gsl_linalg_complex_LU_invert(LU,p,C);
int *pilot=(int)calloc(Nt,sizeof(int));
pilot[k]=1;
bool determine=1;

do
{
for (count1=0;count1<Nt;count1++)
{
gsl_vector_complex_set(C_diag,count1,gsl_matrix_compplex_get(C,count1,count1));
}
int i,k;
gsl_complex Ckk, Cik;
GSL_SET_COMPLEX(&Ckk,0,0);
GSL_SET_COMPLEX(&Cik,0,0);
Ckk=gsl_vector_complex_max(C_diag);
k=gsl_vector_complex_max_index(C_diag);
i=0;
for (count1=0;count1<Nt;count1++)
{
	if(count1==k)
	{
		continue;
	}
	
lamida=-round(gsl_complex_div(gsl_matrix_complex_get(C,count1,k),gsl_matrix_complex_get(C,count1,count1)));
if(lamida==0)
{
	pilot[count1]=1;
}
delta=-gsl_complex_mul_real(gsl_matrix_complex_get(C,count1,count1),gsl_complex_abs2(lamida))-gsl_complex_mul(gsl_complex_conjugate(lamida),gsl_matrix_complex_get(C,count1,k))-gsl_complex_mul(lamida,gsl_matrix_complex_get(C,k,count1));
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
gsl_matrix_complex_get_column(t_k,T,k);
gsl_matrix_complex_get_column(t_i,T,i)
gsl_matrix_complex_get_column(c_k,C,k);
gsl_matrix_complex_get_column(c_i,C,i);
gsl_matrix_complex_get_row(c_kk,T,k);
gsl_matrix_complex_get_row(c_ii,T,i);
gsl_vector_complex_scale(t_i,lamida);
gsl_vector_complex_scale(c_i,lamida)
gsl_vector_complex_scale(c_ii,gsl_complex_conjugate(lamida));
gsl_vector_complex_sub(t_k,t_i);
gsl_vector_complex_sub(c_k,c_i)
gsl_vector_complex_sub(c_kk,c_ii);
gsl_matrix_complex_set_column(t_k,T,k);
gsl_matrix_complex_set_column(c_k,C,k);
gsl_matrix_complex_set_row(c_kk,C,k);
for(count1=0;count1<Nt;count1++)
{
 if(pilot[count1]==0)
 {
	 determine=0;
	 break;
 }
}
}while(determine==0)
//update the matrix pH and T
gsl_matrix_memcpy(LU,T);
gsl_linalg_complex_LU_decomp(LU,p,x);
gsl_linalg_complex_LU_invert(LU,p,T);
gsl_matrix_complex *I=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex_set_identity(I);
gsl_matrix_complex T_update=gsl_matrix_complex_calloc(Nt,Nt);
gsl_matrix_complex pH_update=gsl_matrix_complex_calloc(Nr,Nt);
gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,one,T,I,zero,T_update );
gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,pH,T_update,zero,pH_update);
