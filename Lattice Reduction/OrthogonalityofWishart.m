 Nr=3200;
 Nt=32;
H=complex(normrnd(0,1/sqrt(2),[Nr,Nt]),normrnd(0,1/sqrt(2),[Nr, Nt])); %propagation matrix generation
W=H'*H;    %wishart matrix
W_temp=W;
W_diag=diag(abs(W));
W_min=min(W_diag);
for count1=1:Nr
        W_temp(count1,count1)=0;   
end
W_max=max(max(abs(W_temp)));
W_ratio=W_min/W_max;
sum=1;
for count1=1:Nt
    sum=sum*(norm(H(:,count1),2)^2);
end
det1=abs(det(W));
mesh(1:Nt,1:Nt,abs(W));
orthogonal=1-det1/sum;