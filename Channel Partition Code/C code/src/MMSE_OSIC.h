/*
 * MMSE_OSIC.h
 *
 *  Created on: Dec 9, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "RectangularQAMSlicer.h"
#ifndef MMSE_OSIC_H_
#define MMSE_OSIC_H_
#define LEN sizeof(struct Node)
struct  Node{
	int index;
    struct Node *next;
};
struct Node *create(int Nt){   //create linked list
	struct Node *head, *current,*tail;
	int count;
	for (count=0;count<Nt;count++){
		  current=(struct Node*)calloc(1, LEN);
		  current->index=count;
		  current->next=NULL;
		if(count==0){
          head=current;
          tail=current;
		}
	    tail->next=current;
	    tail=current;
	}
	return head;

}
int get(int k, struct Node *head){   // get the index in the kth struct Node
	int index;
	int count=0;
	struct Node *current=head;
	struct Node *Nulling=NULL;
	if(k==0){
		index=current->index;
		head=current->next;
		free(current);
		current=NULL;
		return index;
	}
	while(count<k-1){
		current=current->next;
		count++;
	}
	Nulling=current->next;
	index=Nulling->index;   //get the index
	current->next=Nulling->next;  //delete Nulling struct Node
	free(Nulling);
	Nulling=NULL;
	return index;
}
//void delete(int k, struct Node *head){   //delete k th struct Node
//	struct Node *current=head;
//	struct Node *Nulling;
//	if(k==0){
//		head=head->next;
//		free(head);
//		return;
//	}
//	int count=0;
//	while(count<k-1){
//		current=current->next;
//		count++;
//	}
//	Nulling=current->next;
//	current->next=Nulling->next;
//	free(Nulling);
//	return;
//}
void MMSE_OSIC(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr,
		double pav, int M, gsl_vector_complex *psymOut){
   int Nr=pH->size1;
   int Nt=pH->size2;
   struct Node *head=create(Nt);
   gsl_matrix_complex *pH_inter=gsl_matrix_complex_calloc(Nr, Nt);
   gsl_matrix_complex_memcpy(pH_inter, pH);
   gsl_vector_complex *preceive_tmp=gsl_vector_complex_calloc(Nr);
   gsl_vector_complex_memcpy(preceive_tmp, preceived);
   gsl_vector_complex_view diag_viewComplex;
   gsl_vector_view diag_viewReal;
   gsl_vector *diag;
   gsl_matrix_complex *pHtemp, *G_pre, *G_preInv, *row_M;
   gsl_permutation *p;
   int *signum=(int*)calloc(1, sizeof(int));
   *signum=1;
   gsl_matrix_complex  *Gmmse=gsl_matrix_complex_calloc(1, Nr);
   gsl_vector_complex  *GmmseR, *colNulling, *colReserve;
   gsl_vector_complex *row;
   GmmseR=gsl_vector_complex_calloc(Nr);
   colNulling=gsl_vector_complex_calloc(Nr);
   colReserve=gsl_vector_complex_calloc(Nr);
   gsl_complex symCurrent;
   gsl_vector_complex *symCurrent_V=gsl_vector_complex_calloc(1);
   gsl_complex alpha, beta1, beta2, beta3;
//   gsl_complex temp;
   int k;
   GSL_SET_COMPLEX(&alpha, 1, 0);
   GSL_SET_COMPLEX(&beta1, pow(snr, -1), 0);
   GSL_SET_COMPLEX(&beta2, 0, 0);
   GSL_SET_COMPLEX(&beta3, -1, 0);
   int count, count1, count2;
   int index;
//   double temp;
   for (count=0;count<Nt;count++){
	   G_pre=gsl_matrix_complex_calloc(Nt-count, Nt-count);
	   G_preInv=gsl_matrix_complex_calloc(Nt-count, Nt-count);
	   gsl_matrix_complex_set_identity(G_pre);
	   p=gsl_permutation_calloc(Nt-count);
	   diag=gsl_vector_calloc(Nt-count);
	   row=gsl_vector_complex_calloc(Nt-count);
	   row_M=gsl_matrix_complex_calloc(1, Nt-count);
//	   if((Nt-count)>1){
//	   gsl_matrix_complex *pHtemp=gsl_matrix_complex_calloc(Nr, Nt-count-1);
//	   }
//	   diag_real=gsl_vector_calloc(Nt-count);
	   gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH_inter,
			   pH_inter, beta1, G_pre);
	   gsl_linalg_complex_LU_decomp(G_pre, p, signum);
	   gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);//calculate the inverse matrix
	   diag_viewComplex=gsl_matrix_complex_diagonal(G_preInv);
	   diag_viewReal=gsl_vector_complex_real(&diag_viewComplex.vector);
	   gsl_vector_memcpy(diag, &diag_viewReal.vector);
	   k=gsl_vector_min_index(diag);   //the current index of the strongest data stream
	   //(with largest post processing SNR
	   gsl_matrix_complex_get_row(row, G_preInv, k);
	   gsl_matrix_complex_set_row(row_M, 0, row);
       if (pH_inter->size1==15){
       	printf("this is the where error happens\n");
       	printf("\n");
       }
	   gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, row_M, pH_inter, beta2, Gmmse);
       gsl_matrix_complex_get_row(GmmseR, Gmmse,0);
       gsl_blas_zdotu(GmmseR, preceive_tmp, &symCurrent);
	   gsl_vector_complex_set(symCurrent_V, 0, symCurrent);  //estimation
	   RectangularQAMSlicer(symCurrent_V, pav, M);
	   //get the orginal index in the list and update the list
//	   struct Node *indexNode;
	   index=get((int)k, head);
//	   index=indexNode->index;
//	   free(indexNode);
	   if (k==0){
		   head=head->next;
	   }
	   gsl_vector_complex_set(psymOut, index, gsl_vector_complex_get(symCurrent_V, 0));
	   //of the chosen symbol
	   //update observation
	   gsl_matrix_complex_get_col(colNulling, pH_inter, k);
	   gsl_vector_complex_scale(colNulling, gsl_vector_complex_get(symCurrent_V, 0));
	   gsl_vector_complex_sub(preceive_tmp, colNulling);
	   //update channel matrix
	   if((Nt-count)==1){
		   gsl_matrix_complex_free(pH_inter);
		   gsl_matrix_complex_free(G_pre);
		   gsl_matrix_complex_free(G_preInv);
		   gsl_permutation_free(p);
		   gsl_vector_free(diag);
		   gsl_vector_complex_free(row);
		   gsl_matrix_complex_free(row_M);



		   break;
	   }

	   pHtemp=gsl_matrix_complex_calloc(Nr, Nt-count-1);
		count2=0;
		for (count1=0;count1<(Nt-count);count1++){
		  if (count1==k){
			  continue;
		  }
		  gsl_matrix_complex_get_col(colReserve, pH_inter, count1);
		  gsl_matrix_complex_set_col(pHtemp, count2, colReserve);
		  count2++;
		  }
		gsl_matrix_complex_free(pH_inter);
		pH_inter=gsl_matrix_complex_calloc(Nr, Nt-count-1);
		gsl_matrix_complex_memcpy(pH_inter, pHtemp);
	    gsl_matrix_complex_free(pHtemp);


	   gsl_matrix_complex_free(G_pre);
	   gsl_matrix_complex_free(G_preInv);
	   gsl_permutation_free(p);
	   gsl_vector_free(diag);
	   gsl_vector_complex_free(row);
	   gsl_matrix_complex_free(row_M);


   }

   gsl_vector_complex_free(preceive_tmp);
   free(signum);
   gsl_vector_complex_free(GmmseR);
   gsl_vector_complex_free(colNulling);
   gsl_vector_complex_free(colReserve);
   gsl_matrix_complex_free(Gmmse);
   gsl_vector_complex_free(symCurrent_V);


	return;

}




#endif /* MMSE_OSIC_H_ */
