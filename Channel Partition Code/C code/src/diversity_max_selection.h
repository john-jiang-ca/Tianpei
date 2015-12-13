/*
 * diversity_max_selection.h
 *
 *  Created on: Dec 11, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "linkedList.h"
#ifndef DIVERSITY_MAX_SELECTION_H_
#define DIVERSITY_MAX_SELECTION_H_
#define LEN sizeof(struct Node)
void diversity_max_selection(gsl_matrix_complex *pH, int N, double snr,
		struct Node *index1_head, struct Node *index2_head,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2){
	int Nr=pH->size1;
	int Nt=pH->size2;
	gsl_combination *subset, *subsetNew;
	gsl_matrix_complex *G_pre, *G_preInv;
	gsl_permutation *p;
	int *signum=(int*)calloc(1, sizeof(int));
	*signum=1;
	gsl_matrix_complex *H2_tmp;
	gsl_vector_complex_view diag_viewComplex;
	gsl_vector_view diag_viewReal;
	gsl_vector *diagMax, *diag;
	gsl_vector_complex *colReserve;
	subset=gsl_combination_calloc(Nt, N);
	int Nu=gsl_sf_fact(Nt)/(gsl_sf_fact(Nt-N)*gsl_sf_fact(N));   //the number of subsets
	int count, count1, count2, count3;
	double max_value;
	p=gsl_permutation_calloc(Nt-N);
	G_pre=gsl_matrix_complex_calloc(Nt-N, Nt-N);
	G_preInv=gsl_matrix_complex_calloc(Nt-N, Nt-N);
	H2_tmp=gsl_matrix_complex_calloc(Nr, Nt-N);
	diagMax=gsl_vector_calloc(Nu);   //store the minimum diagonal value of all the subset
	diag=gsl_vector_calloc(Nt-N);   //the diagonal vector of one subset
	colReserve=gsl_vector_complex_calloc(Nr);
	gsl_complex alpha, beta1, beta2;
	GSL_SET_COMPLEX(&alpha, 1, 0);
	GSL_SET_COMPLEX(&beta1, pow(snr, -1), 0);
	GSL_SET_COMPLEX(&beta2, 0,0);
	int k;
	int pilot=0;
   for (count=0;count<Nu;count++){
	   gsl_matrix_complex_set_identity(G_pre);
      //construct H2_tmp
//	   gsl_matrix_ulong_get_row(subset, subset_M, count);
	   count3=0;
	   for (count1=0;count1<Nt;count1++){
		   pilot=0;
		   for (count2=0;count2<N;count2++){
			   if (count1==gsl_combination_get(subset,count2)){
				   pilot=1;
				   break;
			   }
		   }
		   if(pilot==1){
			   continue;
		   }else{
			   gsl_matrix_complex_get_col(colReserve, pH, count1);
			   gsl_matrix_complex_set_col(H2_tmp, count3, colReserve );
			   count3++;
		   }
	   }
	   gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, H2_tmp, H2_tmp, beta1, G_pre);
	   gsl_linalg_complex_LU_decomp(G_pre, p, signum);
	   gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);
	   diag_viewComplex=gsl_matrix_complex_diagonal(G_preInv);
	   diag_viewReal=gsl_vector_complex_real(&diag_viewComplex.vector);
	   gsl_vector_memcpy(diag, &diag_viewReal.vector);
	   max_value=gsl_vector_max(diag);
	   gsl_vector_set(diagMax, count, max_value);
	   gsl_combination_next(subset);
   }
      k=gsl_vector_min_index(diagMax);
//      gsl_matrix_ulong_get_row(subset, subset_M, k);
     subsetNew=gsl_combination_calloc(Nt, N);
     for (count=0;count<k;count++){  //find the combination number of the best subset
    	 gsl_combination_next(subsetNew);

     }
   //construct pH1, pH2, index1, index2
      index2_head=create(Nt);
      index1_head=split(index2_head, subsetNew);
      count3=0;
      int index;
      for (count1=0;count1<N;count1++){  //construct pH1
    	  index=gsl_combination_get(subsetNew, count1);
    	  gsl_matrix_complex_get_col(colReserve, pH, index);
    	  gsl_matrix_complex_set_col(pH1, count1, colReserve);
      }
      for (count1=0;count1<Nt;count1++){ //construct pH2
     		   pilot=0;
     		   for (count2=0;count2<N;count2++){
     			   if (count1==gsl_combination_get(subsetNew, count2)){
     				   pilot=1;
     				   break;
     			   }
     		   }
     		   if(pilot==1){
     			   continue;
     		   }else{
     			   gsl_matrix_complex_get_col(colReserve, pH, count1);
     			   gsl_matrix_complex_set_col(pH2, count3, colReserve );
     			   count3++;
     		   }
     	   }

  	gsl_combination_free(subset);
  	gsl_combination_free(subsetNew);
  	gsl_matrix_complex_free(G_pre);
  	gsl_matrix_complex_free(G_preInv);
  	gsl_permutation_free(p);
  	free(signum);
  	gsl_matrix_complex_free(H2_tmp);
  	gsl_vector_free(diagMax);
  	gsl_vector_free(diag);
  	gsl_vector_complex_free(colReserve);
	return;
}



#endif /* DIVERSITY_MAX_SELECTION_H_ */
