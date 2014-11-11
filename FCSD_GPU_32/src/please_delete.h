/*
 * please_delete.h
 *
 *  Created on: Sep 29, 2014
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef PLEASE_DELETE_H_
#define PLEASE_DELETE_H_

//elements update
__global__ void linearMge( cuComplex *matrix, int index, int lda) {
//int ty = threadIdx.x;
//int x = blockIdx.x;
//int count1,count2;
//int y = ty + x* blockDim.x;
//cuComplex zero;
//zero.x=0;
//zero.y=0;
////extern __shared__ cuComplex matrix[];
////if(y==0){
////for(count1=0; count1<2*column; count1++){
////	for(count2=0; count2<column; count2++){
////      matrix[IDC2D(count1,count2,lda)]=H[IDC2D(count1,count2,lda)];
////	}
////}
////}
////__shared__ cuComplex multColumn[ LINEAR_BLOCK_SIZE ];
////__shared__ cuComplex matrixPivotValue;
////__shared__ cuComplex matrixRow[ LINEAR_BLOCK_SIZE ];
////__shared__ cuComplex resultPivotValue;
////__shared__ cuComplex resultRow[ LINEAR_BLOCK_SIZE];
////cuComplex newMatrixValue;
//extern __shared__ cuComplex newResultValue[];
////if(y==0){
//// Each block loads the value of the pivot Row to be substracted
//if(x!=index){
//if ( ty == 0 ){
////matrixPivotValue = matrix[ IDC2D( x, index, lda )];
//
//	newResultValue[x]=matrix[ IDC2D( index, x, lda )];
//	matrix[ IDC2D( index, x, lda )]=zero;
//
//}
////multColumn[ ty ] = matrix[ IDC2D( y, index, lda )];
////matrixRow[ ty ] = matrix[ IDC2D( y, x, lda )];
////resultRow[ ty ] = result[ IDC2D( y, x, lda )];
//__syncthreads();
//if ( ty<2*column-index-1 ) {
//matrix[IDC2D(ty+index+1,x,lda)] = complex_sub(matrix[IDC2D(ty+index+1,x,lda)],complex_mulcom(matrix[IDC2D(ty+index+1,index,lda)],newResultValue[x]));
//// Copy to the matrix
//}
//}
//__syncthreads();
////if(ty==0){
////for(count1=0; count1<2*column; count1++){
////	for(count2=0; count2<column; count2++){
////      H[IDC2D(count1,count2,lda)]=matrix[IDC2D(count1,count2,lda)];
////	}
////}
////}
////}
	int ty = threadIdx.x;
	int x = blockIdx.x;
	int y = ty + blockIdx.x * lda;
	 extern __shared__ cuComplex array[ ];
	 cuComplex *multColumn=array;
	 cuComplex *matrixRow=array+2*lda;
//	 cuComplex *multColumn=(cuComplex*)&multColumn1[2*lda*sizeof(cuComplex)];
	 __shared__ cuComplex matrixPivotValue;
//	 cuComplex *matrixRow1=(cuComplex*)array;
//	 cuComplex *matrixRow=(cuComplex*)&matrixRow1[2*lda*sizeof(cuComplex)];

//	__shared__ cuComplex resultPivotValue;
//	__shared__ cuComplex resultRow[ column];
	__shared__ cuComplex newMatrixValue;
	__shared__ cuComplex zero;
	zero.x=0;
	zero.y=0;
	if ( x!=index ) {
	// Each block loads the value of the pivot Row to be substracted
	if ( ty == 0 ){
	matrixPivotValue = matrix[ IDC2D( index, x, lda )];
	printf("the pivot value is:\n");
	printf("%0.4f%+0.4fi ",matrix[ IDC2D( index, x, lda )].x,matrix[ IDC2D( index, x, lda )].y );
//	resultPivotValue = result[ IDC2D( index, x, lda )];
	matrix[ IDC2D( index, x, lda )]=zero;
	}
	__syncthreads();
	multColumn[ ty ] = matrix[ IDC2D( ty, index, lda )];
	matrixRow[ ty ] = matrix[ IDC2D( ty, x, lda )];
//	resultRow[ ty ] = result[ IDC2D( y, x, lda )];
	__syncthreads();
	if ( ty<2*lda-index-1 ) {
//	newMatrixValue =matrix[ IDC2D( ty, x, lda )];
matrix[ IDC2D( ty+index+1, x, lda) ]=complex_sub(matrixRow[ty+index+1],complex_mulcom( multColumn[ty+index+1],matrixPivotValue));
	// Copy to the matrix
//	matrix[ IDC2D( ty, x, lda) ] = newMatrixValue;
	}
	__syncthreads();
	}
}


#endif /* PLEASE_DELETE_H_ */
