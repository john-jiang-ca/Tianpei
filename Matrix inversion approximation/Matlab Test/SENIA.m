function [ W_inverse, W_inter1 ] = SENIA( W )
% This function implements the series expansion (2) Newton Iteration matrix
% inversion appproximation 
%INPUT:
% W: quasi-diagonal matrix
%OUTPUT:
% W_inverse: matrix inversion approximation
n=length(W(1,:));   %length of W
W_diag=zeros(n,1);   %diagonal elements of W   D
W_offdiag=zeros(n); %off diagonal matrix of W  E  
W_diagI=zeros(n,1);  %the elements of the inverse of the diagonal matrix D^(-1)
W_inter1=zeros(n);  %the intermediate matrix D^(-1)E
W_inter2=zeros(n);   %the intermediate matrix W_inter1*D^(-1)
W_inter3=zeros(n);   %the intermediate matrix W_inter2*E
W_inter4=zeros(n);   %the inermediate matrix 
W_diag=diag(W);    %get the diagonal elements of W
W_diagI=1./(W_diag); %get the vector of the elements which is the inverse of the diagonal matrix
W_offdiag=W;
for count=1:n
    W_offdiag(count,count)=0;
end
for count=1:n
    W_inter1(count,:)=W_offdiag(count,:)*W_diagI(count);
end
for count=1:n
    W_inter2(:,count)=W_inter1(:,count)*W_diagI(count);
end
W_inter3=W_inter2*W_offdiag;  
for count=1:n
    W_inter3(count,count)=W_inter3(count,count)+1;
end
W_inter4=-W_inter2;
for count=1:n
W_inter4(count,count)=W_inter4(count,count)+W_diagI(count);
end
W_inverse=W_inter3*W_inter4;
end

