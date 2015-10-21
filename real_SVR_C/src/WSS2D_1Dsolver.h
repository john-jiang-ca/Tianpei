/*
 * WS2D_1Dsolver.h
 *
 *  Created on: Sep 2, 2015
 *      Author: Preston Chen
 */
#include"Public.h"
#ifndef WSS2D_1DSOLVER_H_
#define WSS2D_1DSOLVER_H_
void WSS2D_1Dsolver(gsl_matrix *pH,  //channel propagation matrix
		gsl_vector *symReceived, //received symbol vector
		double SNRd, //signal to noise ratio
		gsl_vector *symbolconstellation, //symbol constellation
		gsl_vector *symOut, //detected symbol vector
		gsl_vector *lambda,//dual variables
		double *Theta, //objective function
		double *G, //duality gap
		int *iteration, //iteration time
		double *MSE //mean square error
		){
	double sigma_tmp, lambda_tmp;
	double lambda_tmp1, lambda_tmp2;
	double delta_tmp1,delta_tmp2;
	gsl_vector *delta=gsl_vector_calloc(Nr);
	gsl_vector *Phi=gsl_vector_calloc(Nr);
	gsl_matrix *K=gsl_matrix_calloc(Nr,Nr);
	gsl_vector *sigma=gsl_vector_calloc(Nr);
//	gsl_vector *sigma_tmp_v=gsl_vector_calloc(Nr);
	double NoiseTerm;
	double best_gain;
	double Theta_tmp;
	int First,Second;
	int count, count1,count2;
	gsl_vector *y_tmp=gsl_vector_calloc(Nr);
	gsl_vector *symOut_tmp=gsl_vector_calloc(Nr);
	gsl_vector *Theta_vector_tmp=gsl_vector_calloc(Nr);
	int sgn[Nr];
	int sgn_1, sgn_2;
//initialization
*iteration=0;
*G=0;
*Theta=0;
*MSE=gsl_blas_dnrm2(symReceived);
gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, pH, pH, 0, K);
gsl_blas_dcopy(symReceived, Phi);
for(count1=0;count1<Nr;count1++){
	if(fabs(gsl_vector_get(Phi,count1))>epsilon){
	*Theta+=pow(fabs(gsl_vector_get(Phi, count1))-epsilon,2);
	}
}
*Theta=-0.5*C*(*Theta);
*G=-2*(*Theta);
*G=*G/fabs(*G+*Theta);
//update lambda
while(*G>tol){
	best_gain=0;
	NoiseTerm=0;
    First=0;
    Second=0;
	for(count=0;count<Nr;count++){
		lambda_tmp1=gsl_vector_get(lambda, count)+(gsl_vector_get(Phi,count)-epsilon*(-1))
				/gsl_matrix_get(K,count, count);
		lambda_tmp2=gsl_vector_get(lambda, count)+(gsl_vector_get(Phi,count)-epsilon*(1))
				/gsl_matrix_get(K,count, count);
//		if(lambda_tmp1<-C){
//				lambda_tmp1=-C;
//			}else if(lambda_tmp1>C){
//				lambda_tmp1=C;
//			}
//		if(lambda_tmp2<-C){
//				lambda_tmp2=-C;
//			}else if(lambda_tmp2>C){
//				lambda_tmp2=C;
//			}
		delta_tmp1=(lambda_tmp1-gsl_vector_get(lambda, count))*(-0.5*(lambda_tmp1
				-gsl_vector_get(lambda,count))*gsl_matrix_get(K,count,count)+gsl_vector_get(Phi,count)
				)-epsilon*(fabs(lambda_tmp1)-fabs(gsl_vector_get(lambda,count)));
		delta_tmp2=(lambda_tmp2-gsl_vector_get(lambda, count))*(-0.5*(lambda_tmp2
				-gsl_vector_get(lambda,count))*gsl_matrix_get(K,count,count)+gsl_vector_get(Phi,count)
				)-epsilon*(fabs(lambda_tmp2)-fabs(gsl_vector_get(lambda,count)));
		if(delta_tmp1>delta_tmp2){
			lambda_tmp=lambda_tmp1;
			sgn[count]=-1;
			gsl_vector_set(delta,count ,delta_tmp1);
		}else{
			lambda_tmp=lambda_tmp2;
			sgn[count]=1;
			gsl_vector_set(delta, count, delta_tmp2);
		}


//		sigma_tmp=lambda_tmp-gsl_vector_get(lambda,count);
//		gsl_vector_set(sigma, count, sigma_tmp);


	}
#ifdef DEBUG
	printf("the gain of objective function are\n");
	for(count2=0;count2<Nr;count2++){
		printf("%f ", gsl_vector_get(delta,count2));
	}
	printf("\n");
#endif
//sorting
//find the first coordinate that can generate the first and
//the second largest objective function gain
	First=0;
	best_gain=gsl_vector_get(delta,0);
	for(count2=1;count2<Nr;count2++){
		if(gsl_vector_get(delta, count2)>best_gain){
			best_gain=gsl_vector_get(delta,count2);
			First=count2;
		}
	}

//find the second coordinate
	Second=-1;
	best_gain=0;
	for(count2=0;count2<Nr;count2++){
		if(count2==First){
			continue;
		}
		if(gsl_vector_get(delta, count2)>best_gain){
			best_gain=gsl_vector_get(delta,count2);
			Second=count2;
		}
	}
//update lambda, sigma, Phi, Theta and duality gap
//	lambda_tmp=gsl_vector_get(lambda,First)+gsl_vector_get(sigma, First);
//	gsl_vector_set(lambda, First, lambda_tmp);
//	lambda_tmp=gsl_vector_get(lambda,Second)+gsl_vector_get(sigma, Second);
//	gsl_vector_set(lambda, Second, lambda_tmp);

//call 2-D solver
//update lambda and sigma
	sgn_1=sgn[First];
	sgn_2=sgn[Second];
	lambda_tmp=gsl_vector_get(lambda, First)+(gsl_vector_get(Phi,First)
			*gsl_matrix_get(K,Second,Second)-gsl_vector_get(Phi, Second)
			*gsl_matrix_get(K,First,Second)-epsilon*(sgn_1*gsl_matrix_get(K,Second,Second)
	-sgn_2*gsl_matrix_get(K,First, Second)))/(gsl_matrix_get(K,First,First)
					*gsl_matrix_get(K,Second,Second)-pow(gsl_matrix_get(K,First,Second),2));
	sigma_tmp=lambda_tmp-gsl_vector_get(lambda, First);
	gsl_vector_set(sigma,First,sigma_tmp);
	gsl_vector_set(lambda, First, lambda_tmp);
	lambda_tmp=gsl_vector_get(lambda, Second)+(gsl_vector_get(Phi,Second)
			*gsl_matrix_get(K,First,First)-gsl_vector_get(Phi, First)
			*gsl_matrix_get(K,First,Second)-epsilon*(sgn_2*gsl_matrix_get(K,First,First)
	-sgn_1*gsl_matrix_get(K,First, Second)))/(gsl_matrix_get(K,First,First)
					*gsl_matrix_get(K,Second,Second)-pow(gsl_matrix_get(K,First,Second),2));
	gsl_vector_get(lambda, Second);
	sigma_tmp=lambda_tmp-gsl_vector_get(lambda, Second);
	gsl_vector_set(sigma,Second,sigma_tmp);
	gsl_vector_set(lambda, Second, lambda_tmp);

//update Phi
	double Phi_tmp;
	for(count2=0;count2<Nr;count2++){
		Phi_tmp=gsl_vector_get(Phi,count2)
				-gsl_vector_get(sigma,First)*gsl_matrix_get(K, count2, First)
				-gsl_vector_get(sigma,Second)*gsl_matrix_get(K, count2, Second);
		gsl_vector_set(Phi,count2, Phi_tmp);
	}
#ifdef DEBUG
	printf("\n");
	printf("the Phi after updating are:\n");
	for(count2=0;count2<Nr;count2++){
		printf("%f ", gsl_vector_get(Phi,count2));
	}
	printf("\n");
#endif

//	Theta_tmp=-(1/2)*lambda'*K*lambda+y'*lambda-epsilon*norm(lambda, 1);
//	Theta(iteration)=Theta_tmp-0.5*NoiseTerm;
	for(count1=0;count1<Nr;count1++){
		if(fabs(gsl_vector_get(Phi, count1))-epsilon>0){
			NoiseTerm+=pow((double)(fabs(gsl_vector_get(Phi, count1))-epsilon),2);
		}
	}

	double *Theta_tmp1=(double*)malloc(sizeof(double));
	double *Theta_tmp2=(double*)malloc(sizeof(double));
	double Theta_tmp3=0;
	gsl_blas_dgemv(CblasNoTrans, -0.5, K, lambda, 0, Theta_vector_tmp);
	gsl_blas_ddot(lambda, Theta_vector_tmp, Theta_tmp1);
	gsl_blas_ddot(symReceived, lambda, Theta_tmp2);
	Theta_tmp3=epsilon*gsl_blas_dasum(lambda);
	Theta_tmp=*Theta_tmp1+*Theta_tmp2-Theta_tmp3-0.5*C*NoiseTerm;
	 (*iteration)++;
	*Theta=Theta_tmp;
	double G_tmp=*Theta_tmp2-Theta_tmp3-2*Theta_tmp;
//	printf("G_tmp now is %f'n ", G_tmp);
	*G=G_tmp/fabs(G_tmp+*Theta);
#ifdef DEBUG
	printf("G now is %f\n", *G);
	printf("Theta now is %f\n", *Theta);
#endif
    free(Theta_tmp1);
    free(Theta_tmp2);

}

//calculate detected symbol
gsl_blas_dgemv(CblasTrans, 1, pH, lambda, 0, symOut);
//rounding
double d=sqrt((double)3/(Nt*(pow(M,2)-1)));
for(count1=0;count1<Nt;count1++){
	if(M==4){
		if(gsl_vector_get(symOut, count1)<=-2*d){
         gsl_vector_set(symOut,count1, -3*d);
		}else if (gsl_vector_get(symOut ,count1)>-2*d&&gsl_vector_get(symOut, count1)<=0){
 gsl_vector_set(symOut, count1, -d);
		}else if(gsl_vector_get(symOut,count1)>0&&gsl_vector_get(symOut, count1)<=2*d){
gsl_vector_set(symOut, count1, d);
		}else if(gsl_vector_get(symOut ,count1)>2*d){
gsl_vector_set(symOut, count1 , 3*d);
		}

	}else if(M==8){
		if(gsl_vector_get(symOut, count1)<=-6*d){
         gsl_vector_set(symOut,count1, -7*d);
		}else if (gsl_vector_get(symOut ,count1)>-6*d&&gsl_vector_get(symOut, count1)<=-4*d){
 gsl_vector_set(symOut, count1, -5*d);
		}else if(gsl_vector_get(symOut,count1)>-4*d&&gsl_vector_get(symOut, count1)<=-2*d){
gsl_vector_set(symOut, count1, -3*d);
		}else if(gsl_vector_get(symOut ,count1)>-2*d&&gsl_vector_get(symOut, count1)<=0){
gsl_vector_set(symOut, count1 , -d);
		}else if(gsl_vector_get(symOut, count1)>6*d){
	         gsl_vector_set(symOut,count1, 7*d);
			}else if (gsl_vector_get(symOut ,count1)>4*d&&gsl_vector_get(symOut, count1)<=6*d){
	 gsl_vector_set(symOut, count1, 5*d);
			}else if(gsl_vector_get(symOut,count1)>2*d&&gsl_vector_get(symOut, count1)<=4*d){
	gsl_vector_set(symOut, count1, 3*d);
			}else if(gsl_vector_get(symOut ,count1)>0&&gsl_vector_get(symOut, count1)<=2*d){
	gsl_vector_set(symOut, count1 , d);
			}

	}
}

//calculate MSE
gsl_blas_dcopy(symReceived, y_tmp);
gsl_blas_dgemv(CblasNoTrans, 1 , pH, symOut, 0, symOut_tmp);
gsl_vector_sub(y_tmp,symOut_tmp);
*MSE=gsl_blas_dnrm2(y_tmp);

gsl_vector_free(delta);
gsl_vector_free(Phi);
gsl_matrix_free(K);
gsl_vector_free(sigma);
gsl_vector_free(y_tmp);
gsl_vector_free(symOut_tmp);
gsl_vector_free(Theta_vector_tmp);
//gsl_vector_free(sigma_tmp_v);
return;
}




#endif /* WSS1D_2DSOLVER_H_ */
