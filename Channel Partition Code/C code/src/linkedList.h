/*
 * linkedList.h
 *
 *  Created on: Dec 11, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef LINKEDLIST_H_
#define LINKEDLIST_H_

#define LEN sizeof(struct Node)
struct Node{
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


struct Node *split(struct Node *index2_head, gsl_combination *subset){
	//split the complete linked list index1 into two sub linked list index1 and index2
	struct Node *index1_head;
	struct Node *current, *first, *Nulling, *index1_tail;
	current=create(1);
	first=current;
	current->next=index2_head;
	int count;
	while(current->next!=NULL){
		if(current->next->index==gsl_combination_get(subset, count)){
			Nulling=current->next;
			current->next=Nulling->next;
			Nulling->next=NULL;
			if(count==0){
				index1_head=Nulling;
				index1_tail=Nulling;
			}else{
			index1_tail->next=Nulling;
			index1_tail=index1_tail->next;
			}
			count++;
		}
		current=current->next;
	}
	free(first);
	return index1_head;

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



#endif /* LINKEDLIST_H_ */
