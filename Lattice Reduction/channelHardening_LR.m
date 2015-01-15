function [ k,W_ratio,W_LR_ratio,orthogonal,orthogonal_LR,Var,Var_o ] = channelHardening_LR( Nr,Nt )
%% Test the effect of lattice reduction to the channel hardening process
%INPUT
%Nr: number of receive antennas
%Nt: number of transmit antennas
%OUTPUT
%k: feedback
%the 3-D figure shows the results

%Nr=8;
%Nt=8;
H=complex(normrnd(0,1/sqrt(2),[Nr,Nt]),normrnd(0,1/sqrt(2),[Nr, Nt])); %propagation matrix generation
H_temp=H;
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
orthogonal=1-det(W)/sum;
C=eye(Nt,Nt)/W;
%% COMPLEX LATTICE REDUCTION ALGORITHM (CLLL)
InnerVec=zeros(1,Nt);   %the vector used to store the inner products
U=eye(Nt,Nt);     %Unitary matrix
Uindex=zeros(Nr,Nt); %the index of Unitary matrix
for count1=1:Nt
    hj=H(:,count1);
    InnerVec(count1)=hj'*hj;
end
for count1=1:Nt
for count2=count1+1:Nt
    a=0;
    for count3=1:count1-1
        a=a+Uindex(count1,count3)'*Uindex(count2,count3)*InnerVec(count3);
    end
    Uindex(count2,count1)=(H(:,count2)'*H(:,count1)-a)/InnerVec(count1);
    InnerVec(count2)=InnerVec(count2)-(norm(Uindex(count2,count1),2)^2)*InnerVec(count1);
end
end

count4=2;
while (count4<=Nt)
    if(real(Uindex(count4,count4-1))>1/2||imag(Uindex(count4,count4-1)>1/2))
        c=round(Uindex(count4,count4-1));
        H(:,count4)=H(:,count4)-c*H(:,count4-1);
        U(:,count4)=U(:,count4)-c*U(:,count4-1);
        for count5=1:count4-1
            Uindex(count4,count5)=Uindex(count4,count5)-c*Uindex(count4-1,count5);
        end
    end
    if InnerVec(count4)<(1/4-norm(Uindex(count4,count4-1),2)^2)*InnerVec(count4-1); %check swaping condition
       
       
        %update acccording to formulas in paper
       % v1=H(:,count4);
       % v2=H(:,count4-1);
        f=Uindex(count4,count4-1);
        Uindex(count4,count4-1)=Uindex(count4,count4-1)*InnerVec(count4-1)/(InnerVec(count4)+Uindex(count4,count4-1)*Uindex(count4,count4-1)'...
            *InnerVec(count4-1));
        ff=zeros(Nt);
        for count=count4+1:Nt
            ff(count)=Uindex(count,count4-1);
            Uindex(count,count4-1)=Uindex(count,count4-1)*Uindex(count4,count4-1)-Uindex(count,count4)*InnerVec(count4)...
                /(InnerVec(count4)+f*f'*InnerVec(count4-1));
             Uindex(count,count4)=ff(count)-Uindex(count,count4)*f;
           
        end
       % ff(count4)=f;
       % for count=count4:Nt 
        %end
        H(:,[count4-1;count4])=H(:,[count4;count4-1]);
        Uindex([count4-1;count4],1:count4-2)=Uindex([count4;count4-1],1:count4-2);
        U(:,[count4-1,count4])=U(:,[count4,count4-1]);
        count4=max(2,count4-1);
    else
        for count5=count4-2:-1:1
            if (real(Uindex(count4,count5))>1/2||imag(Uindex(count4,count5)>1/2))
           c=round(Uindex(count4,count5));
        H(:,count4)=H(:,count4)-c*H(:,count5);
        U(:,count4)=U(:,count4)-c*U(:,count5);
        for count6=1:count5
            Uindex(count4,count6)=Uindex(count4,count6)-c*Uindex(count5,count6);
        end
            end
        end
        count4=count4+1;
            
    end
end
    W_LR=H'*H;
    W_LR_temp=W_LR;
W_LR_diag=diag(abs(W_LR));
W_LR_min=min(W_LR_diag);
for count1=1:Nr
        W_LR_temp(count1,count1)=0;   
end
W_LR_max=max(max(abs(W_temp)));
W_LR_ratio=W_LR_min/W_LR_max;
sum=1;
for count1=1:Nt
    sum=sum*(norm(H(:,count1),2)^2);
end
orthogonal_LR=1-det(W_LR)/sum;
C_o=eye(Nt,Nt)/W_LR;
Var=var(diag(C));
Var_o=var(diag(C_o));
%% PLOT FIGURES
%x=1:Nt;
%y=1:Nt;
%figure (1)
%mesh(x,y,abs(W_LR)),title('LR based wishart matrix'),xlabel('receive antennas index'),ylabel('transmit antennas index');  %plot the diagonal 3-D figure of wishart matrix after LR
%figure (2)
%mesh(x,y,abs(W)),title('original wishart matrix'),xlabel('receive antennas index'),ylabel('transmit antennas index');  %plot the diagonal 3-D figure of original wishart matrix
%figure (3)
%mesh(1:Nt,1:Nr,abs(H_temp-H)),title('The change between original matrix and lattice reducted matrix');
%figure (4)
%mesh(x,y,abs(C)),title('PEP original mesh');
%figure (5)
%mesh(x,y,abs(C_o)),title('PEP orthogonalization mesh');
%figure (6)
%plot(x,diag(C)),title('PEP original plot');
%figure (7)
%plot(x,diag(C_o)),title('PEP orthogonalization plot');
if (abs(orthogonal)-abs(orthogonal_LR)>=0)
k=1;
else
k=0;
end
%k=1;
end

%%  Using for testing
% for count=1:1e4
%[ k,W_ratio,W_LR_ratio,orthogonal,orthogonal_LR,Var,Var_o ] = channelHardening_LR( 4,4 );
%o_v(count)=orthogonal;
%o_v_lr(count)=orthogonal_LR;
%Var_v(count)=Var;
%Var_o(count)=Var_o;
%end
%o_v=round(o_v*1e2)/1e2;
%o_v_lr=round(o_v_lr*1e2)/1e2;
%table1=tabulate(o_v(:));
%table2=tabulate(o_v_lr(:));
%figure (1), plot(table1(:,1),table1(:,3)),hold on;
%plot(table2(:,1),table2(:,3),'-r'),hold off;

