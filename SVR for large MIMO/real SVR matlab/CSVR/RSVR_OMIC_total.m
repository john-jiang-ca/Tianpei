function [ symOut ] = RSVR_OMIC_total(y, H, symOut_prev, SNRd, M, pav, stage, maxStage,W,G, symConstell, order,C, tol, epsilon)
%Multistage Interference Cancellation (MIC) recursive algorithm
%   Detailed explanation goes here
% W is weight matrix with the diagonal elements represent the weight of
% each estimation
%is W is an identity matrix, that this is the original version MIC
%the weakest data stream is detected firsly
Nr=length(H(:,1)); %number of receive antennas
Nt=length(H(1,:));  %number of transmit antennas
I=eye(Nt);
% maxStage=4;
%% First stage estimation by MMSE
% if(Nstage==0)
symOut_current=zeros(Nt,1);
if(stage==1)
INV=(I/(H'*H+SNRd^(-1)*I));
G=INV*H';
D=diag(INV);
[value, order]=sort(D,'ascend');
% H_tmp=H;
% G_tmp=INV;
% order=zeros(Nt,1);
% list=1:Nt;
% for count=1:Nt-1
%     [value, i]=min(diag(G_tmp));
%     order(count)=list(i);
%     H_tmp(:,i)=[];
%     list(i)=[];
%     I=eye(Nt-count);
%     G_tmp=(I/(H_tmp'*H_tmp+SNRd^(-1)*I));
% end
% order(Nt)=list(1);
% symOut_current=G*y;
% [symOut_current]=Rectangular_QAM_slicer(symOut_current, M, pav);
H_r=[real(H), -imag(H); imag(H), real(H)];
y_r=[real(y);imag(y)];
[ symOut_current] = RSVR( H_r,  y_r, SNRd,  M, pav,C,tol, epsilon );
symOut_prev=symOut_current;
end


%    symOut=symI;
%    return;
% Nstage=Nstage+1;
% end
% if(Nstage==MAXstage)
%     return;
% end
%% Multistage-stage PIC
% symOut_current=zeros(Nt,1);
%% weighted PIC
% [W]=weightCal(G,SNRd, symOut_prev, symConstell);
W=eye(Nt,Nt);
symOut_tmp=symOut_prev;
for i=1:Nt
    count=order(i);
    W_t=W;
    W_t(count,count)=0;
    y_tmp=y-H*W_t*symOut_tmp;
%     symOut_tmp(count)=((H(:,count)')/(H(:,count)'*H(:,count)+SNRd^(-1)))*y_tmp;
y_tmp_r=[real(y_tmp); imag(y_tmp)];
H_tmp_r=[real(H(:,count)), -imag(H(:,count)); imag(H(:,count)), real(H(:,count))];
     [symOut_tmp(count)]=RSVR_sub(H_tmp_r, y_tmp_r, SNRd, M, pav, C, tol, epsilon);
%     [symOut_tmp(count)]=Rectangular_QAM_slicer(symOut_tmp(count),M, pav);
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
[symOut]=RSVR_OMIC_total(y, H, symOut_current, SNRd, M, pav, stage+1, maxStage, W,G, symConstell,order, C,tol, epsilon);
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