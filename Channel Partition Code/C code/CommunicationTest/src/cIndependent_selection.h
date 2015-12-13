/*
 * cIndependent_selection.h
 *
 *  Created on: Dec 12, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "linkedList.h"
#ifndef CINDEPENDENT_SELECTION_H_
#define CINDEPENDENT_SELECTION_H_
void cIndependent_selection(gsl_matrix_complex *pH, int N, gsl_rng *pr,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2,
		struct Node *index1_head, struct Node *index2_head){
	int Nr=pH->size1;
	int Nt=pH->size2;
	int Nu=gsl_sf_fact(Nt)/(gsl_sf_fact(Nt-N)*gsl_sf_fact(N));   //the number of subsets
    gsl_combination *subset;
    gsl_vector_complex *colReserve;
    colReserve=gsl_vector_complex_calloc(Nr);
    int k=gsl_rng_uniform_int(pr, Nu);   //randomly choose a combination subset
    int index;
    int count=0;
    int count1,count2;
    int pilot;
    subset=gsl_combination_calloc(Nt, N);
    while (count<k){
    	gsl_combination_next(subset);
    	count++;
    }
    index2_head=create(Nt);
    index1_head=split(index2_head, subset);
    for(count=0;count<N;count++){
    	index=gsl_combination_get(subset, count);
    	gsl_matrix_complex_get_col(colReserve, pH, index);
    	gsl_matrix_complex_set_col(pH1, count, colReserve);
    }
    count2=0;
    for (count=0;count<Nt;count++){
    	pilot=0;
    	for(count1=0;count1<N;count1++){
    		if(count==gsl_combination_get(subset, count1)){
    			pilot=1;
    			break;
    		}
    	}
    	if(pilot==0){
    	 gsl_matrix_complex_get_col(colReserve, pH, count);
    	 gsl_matrix_complex_set_col(pH2, count2, colReserve);
    	 count2++;
    	}
    }
    gsl_combination_free(subset);
    gsl_vector_complex_free(colReserve);
	return;
}




#endif /* CINDEPENDENT_SELECTION_H_ */
