/*
 ============================================================================
 Name        : CC_practice.c
 Author      : Tianpei Chen
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <>
int fix(float a){
	int answer;
	answer=(int)(a<0? (a+0.5): (a-0.5));
	return answer;

}


void reverse1(char str[]){
printf("the original string is %s\n", str);
int count; //counter
size_t length;
length=strlen(str);
char charTemp;
int end=(int)((length-1)/2);
printf("%d\n",end);
for(count=0;count<=end;count++){
	charTemp=str[count];
	int first=count;
	int last=(length-1)-count;
	charTemp=str[first];
	str[first]=str[last];
	str[last]=charTemp;
}
printf("the reverse string is %s\n", str);
}

void reverse2(char *str){
	printf("the original string is %s\n", str);
	char *end, *strtmp, temp;
	end=str;
	strtmp=str;    //used to store the string pointer
	while(*end){ //find the end of the string
		++end;
	}
	--end;  //go back for one bit  because the last bit is NULL
	/*swap the string elements*/
	while(end>strtmp){
		temp=*strtmp;
		*strtmp++=*end;
		*end--=temp;
	}
	// the pointer string is still point to the begining of the string
	printf("the reverse string is %s\n", str);
}



int main(){
	int fix(float a);
	void reverse1(char str[]);
	void reverse2(char *str);
	char str1[8]="tianpei";
	char *str2=str1;
	reverse1(str1);
	reverse2(str2);
	return 0;
}
