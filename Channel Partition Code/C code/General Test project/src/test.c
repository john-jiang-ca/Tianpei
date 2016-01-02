/*
 ============================================================================
 Name        : test.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(void) {
	double a= 3.3;
	double b;
	b=ceil(a);
    int Nr=16;
    int Nt=20;
    double Ntmp=sqrt(20);
//    double Ntmp1=sqrt(Nr+(double)pow((Nr-Nt),2.0)/(double)4);
//    int Ntmp1=(double)(Nr-Nt)*0.5;
    double Ntmp1=0.25*pow((Nr-Nt),2.0);
	int N=ceil(sqrt((double)Nr+(double)pow((Nr-Nt),2.0)/(double)4)-0.5*(double)(Nr-Nt)-1);
	printf("the result for ceil is %f\n", b);
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */
	return EXIT_SUCCESS;
}
