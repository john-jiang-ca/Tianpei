/*
 * SVR_DETECTOR.h
 *
 *  Created on: Jul 15, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SVR_DETECTOR_H_
#define SVR_DETECTOR_H_
#include "common.h"
/*
 *  Parameters configuration
 */
//#define  C  1e-3;      //Penalize weight for noise
//#define epsilon 1e-3; //Training precision
typedef struct{
	gsl_complex data;
	int index;

}Nalpha_Node;
/*
 * inner function
 */
//int TestData(gsl_complex alpha, int index, gsl_vector_complex *Error, gsl_vector_complex *eta);
//int takestep(gsl_complex alpha, int index, gsl_vector_complex *Error, gsl_vector_complex *eta);
//Nalpha_Node checkNonbound(gsl_vector_complex *alpha);

int Initialization(gsl_vector *alpha,//Lagrange multiplier
		gsl_vector *beta,//Lagrange multiplier
		const gsl_vector_complex *preceived,//received symbol vector(output of trainning data set)
		gsl_vector_complex *phi, // update parameter
		double S_C_real, //real stopping parameter
		double S_C_imag,  //imaginary stopping parameter
		int method //label which start strategy is used
		);
void 	WSS2_1Dsolver(gsl_vector *l_m,  //Lagrange Multiplier vector
        gsl_vector_complex *phi,  //update parameter
        gsl_matrix *R_Kernel,    //real kernel matrix
        int F_i,  //index of first maximum Lagrange multiplier
        int S_i, //index of second maximum Lagrange multiplier
		int model //determine this routine is for real part (0) or imaginary part (1)
);
void SVR_DETECTOR(const gsl_vector_complex *preceived, const gsl_matrix_complex *pH, float SNRb, gsl_vector_complex *psymout, int start_M,int selection_M){
//gsl_vector  *eta=gsl_vector_calloc(Nt); //the diagonal components of Wishart matrix
//gsl_vector_complex *Error=gsl_vector_complex_calloc(Nr);  //the error e=y-hx
int count1, count2, count3; //counters
double S_pilot=0; //stoping criteria parameter
gsl_vector *alpha=gsl_vector_calloc(2*Nr);
//gsl_vector *alpha_hat=gsl_vector_calloc(Nr);
gsl_vector *beta=gsl_vector_calloc(2*Nr);
//gsl_vector *beta_hat=gsl_vector_calloc(Nr);  // Lagrange multipliers
gsl_matrix_complex *K_complex=gsl_matrix_complex_calloc(Nr,Nr);
gsl_blas_zherk(CblasUpper,CblasNoTrans, 1, pH,0, K_complex);  //calculate complex kernel matrix
gsl_matrix *K=gsl_matrix_calloc(Nr,Nr);   //kernel matrix
for(count1=0;count1<Nr;count1++){
	for(count2=0;count2<=count1;count2++){
		gsl_matrix_set(K,count1,count2,gsl_matrix_complex_get(K_complex,count1,count2).dat[0]);
	}
}   //get the real kernel
//const gsl_matrix_complex_view K=gsl_matrix_complex_real(K_complex);

gsl_vector_complex *phi=gsl_vector_complex_calloc(Nr); //phi matrix (update parameter)
gsl_vector *delta_psi=gsl_vector_calloc(Nr); //gain of objective function
double S_C_real=0; //stopping criteria(real)
double S_C_imag=0; //stopping criteria(imaginary)
// initialization
int M=Constellationsize;  //the size of modulation constellation
Initialization(alpha,beta,preceived,phi,S_C_real,S_C_imag, start_M);
double L=0;
double H=C*double(1)/double(SNRb*log2(M));
//alpha is chosen by the random value in the interval [0, C sigma], where sigma is the standard
//derivative of the noise;
int count1, count2, count3; //counters
int i, j;
count2=0; //iteration times of real part
int model_real=0;
int model_imag=1;
while(stopping criteria (real) not satisfied){
	WSS2_1Dsolver(alpha,phi,K,i,j,0);
	TakeStep(alpha,phi,K,0,i,j);
	count2++;

}
i=j=0;
count3=0; //iteration times of imaginary part
while(stopping criteria(imaginary) not satisfied){
	WSS2_1Dsolver(beta, phi, K , i, j,1);
	TakeStep(alpha,phi,K,1,i,j);
	count3++;
}

//gsl_complex e_tmp,result_tmp;
// gsl_vector_complex *hi=gsl_vector_complex_calloc(Nt);
// gsl_vector_complex *hj=gsl_vector_complex_calloc(Nt);
//for(count1=0;count1<Nr;count1++){
//	GSL_SET_COMPLEX (&e_tmp, 0, 0);
//	GSL_SET_COMPLEX(&result_tmp,0,0);
//	gsl_matrix_complex_get_row(hi,pH,count1);
//for(count2=0;count2<Nr;count2++){
//   gsl_matrix_complex_get_row(hj,pH,count2);
//	gsl_blas_cdotu( hi, hj,result_tmp);
//	e_tmp+=(gsl_vector_complex_get(alpha,count2)-gsl_vector_complex_get(alpha,count2+Nr))*result_tmp;
//	gsl_vector_complex_set(Error, gsl_vector_complex_get(preceive,count1)-e_tmp, count1);
//}
//}
//gsl_vector_complex_free(hi);
//gsl_vector_complex_free(hj);


/*
 * main routine
 */
//void outterloop(alpha, eta, Error, pH, psymout){
////	begine procedure
//Nalpha_Node Nalpha[Nr];
//Nalpha_Node=checkNonbound(alpha);  //find the non bound alpha and store them in Nalpha
//int Numchanged=0; // the number of alpha changed
//int examineAll=1;  // whether check the whole alpha set all just the non bound alpha
//while (Numchanged>0||examineAll==0){
//	Numchanged=0;
//	if(examineAll){
//		for(count1=0;count1<2*Nr;count1++){
//			gsl_complex tmp;
//			tmp=gsl_vector_complex_get(alpha(count1));
//         Numchanged+=TestData(tmp,count1, Error, eta);
//		}
//	}else{
//		for(count1=0;count1<;count1++){
//			gsl_complex tmp;
//			tmp=Nalpha[count1].data;
//			Numchanged+=TestData(tmp, Nalpha[count1].index, Error, eta);
//
//		}
//	}
//	if(examineAll==1){
//		examineAll=0;
//	}else if(Numchanged==0){
//		examineAll=1;
//	}
//}
//
//}






}





#endif /* SVR_DETECTOR_H_ */
