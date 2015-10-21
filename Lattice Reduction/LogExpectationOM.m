%This programme visualize the logarithmic expectation of orthogonality
%measure
%Preston Chen
%Feb/24/2015
E=zeros(40,20);

for Nt=1:5:100
    for Nr=Nt:5:200
        E((Nr-1)/5+1,(Nt-1)/5+1)=Nt*log(2)/log(exp(1))-Nt*psi(Nr);
        for i=1:Nt
        E((Nr-1)/5+1,(Nt-1)/5+1)=E((Nr-1)/5+1,(Nt-1)/5+1)+psi((Nr+1-i)/2);
        end
    end
end
figure (1)
mesh(1:5:100,1:5:200,E),xlabel('Nt'),ylabel('Nr'),zlabel('Orthogonality Measurement');