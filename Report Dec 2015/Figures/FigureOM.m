%% this routine calcultate the theoretical logarithmic expectation of Orthogonality measure (OM), as well as the 
% the empirical estimation of the logarithmic expecation of OM.
% the channel matrix sizes considered here is 5<=Nr<=200, and 5<=Nt<=Nr
count1=0;
count2=0;
count3=0;
Crealization=1e3;
e=exp(1);
Size=20;
Et=zeros(Size);
Eem=zeros(Size);
count=0;
for Nr=5:5:5*Size
    count1=count1+1;
    count2=0;
    for Nt=5:5:Nr
        count2=count2+1;
        Et(count1, count2)=0;
        for i=1:Nt
        Et(count1,count2)=Et(count1, count2)+(psi(Nr-i+1)-psi(Nr));% calculate the theoretical logarithmic expectation of OM
        end
         Eem(count1, count2)=0;
         
        for count3=1:Crealization   %calculate the empirical estimation of the logarithmic expectation of OM 
           H=complex(normrnd(0,1/sqrt(2),[Nr,Nt]),normrnd(0,1/sqrt(2),[Nr, Nt]));  %generate the channel matrix
           hNorm=1;
           for count4=1:Nt
               hNorm=hNorm*norm(H(:,count4))^2;
           end
           Eem(count1, count2)=Eem(count1, count2)+log(det(H'*H)/(hNorm))/log(e);
        end
        Eem(count1, count2)=Eem(count1,count2)/Crealization;    %taking average of the 
        count=count+1;
    end
end

V_M=Eem-Et;
V=0;    %% the variance between E(ln(OM))t and E(ln(OM))em
for i=1:Size
    for j=1:Size
        V=V+V_M(i,j)^(2);
    end
end
V=V/(count);
figure (1)
mesh(5:5:5*Size,5:5:5*Size,Et),xlabel('Nt'),ylabel('Nr'),zlabel('Orthogonality Measurement(t)');
figure (2)
mesh(5:5:5*Size,5:5:5*Size,real(Eem)),xlabel('Nt'),ylabel('Nr'),zlabel('Orthogonality Measurement(em)');