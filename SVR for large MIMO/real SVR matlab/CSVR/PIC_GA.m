function [ symOut ] = PIC_GA(y, H, symOut_prev, SNRd, M, pav, stage, maxStage, symConstell, order, W)
%Multistage Interference Cancellation (MIC) recursive algorithm
%   Detailed explanation goes here
% W is weight matrix with the diagonal elements represent the weight of
% each estimation
%is W is an identity matrix, that this is the original version MIC
%the weakest data stream is detected firsly
Nr=length(H(:,1)); %number of receive antennas
Nt=length(H(1,:));  %number of transmit antennas
% if(stage==1)
% % I=eye(Nt);
% W=eye(Nt,Nt);
% G=H'*H;
% [value, order]=sort(diag(G),'ascend');
% end
symOut_tmp=symOut_prev;
E=zeros(M,1);

for i=1:Nt
    count=order(i);
    if(i==1)
     W_t=W;
    W_t(count,count)=0;
    y_tmp=y-H*W_t*symOut_tmp;
    else
    y_tmp=E_matrix(:,index1)+H(:,count)*symOut_tmp(count);
    end
 E_matrix=y_tmp*ones(1,M)-H(:,count)*symConstell.';
    for count1=1:M
        E(count1)=norm(E_matrix(:,count1));
    end
    [value1, index1]=min(E);
    symOut_tmp(count)=symConstell(index1);
 
%     symI(i)=sym2(i);
end
symOut_current=symOut_tmp;
% [symOut_current]=Rectangular_QAM_slicer(symOut_current,M, pav);
%   symOut=sym2;
%   return;
if(stage==maxStage)
    symOut=symOut_current;
    return
end
[symOut]=PIC_GA(y, H, symOut_current, SNRd, M, pav, stage+1, maxStage, symConstell, order,W);
end


%% adaptive  weight function 
function [W]=weightCal(G,SNRd, symOut_prev, symConstell)
% d=sqrt(3*pav/(2*(M-1)));
Nt=length(G);
noiseV=sqrt(1/SNRd);
W=zeros(Nt,Nt);
for count1=1:Nt
   if(real(symOut_prev(count1))<0)
       p1=qfunc(real(symOut_prev(count1))/(noiseV*norm(G(count1,:))));
   else
       p1=qfunc(-real(symOut_prev(count1))/(noiseV*norm(G(count1,:))));
   end
   
   if(imag(symOut_prev(count1))<0)
       p2=qfunc(imag(symOut_prev(count1))/(noiseV*norm(G(count1,:))));
   else
       p2=qfunc(-imag(symOut_prev(count1))/(noiseV*norm(G(count1,:))));
   end
   
     W(count1,count1)=p1*p2;
end
end