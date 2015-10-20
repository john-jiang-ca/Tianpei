function [ symOut ] = MIC_Recursive(y, H, symOut_prev, SNRd, M, pav, stage, maxStage,W,G, symConstell)
%Multistage Interference Cancellation (MIC) recursive algorithm
%   Detailed explanation goes here
% W is weight matrix with the diagonal elements represent the weight of
% each estimation
%is W is an identity matrix, that this is the original version MIC
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
symOut_current=G*y;
[symOut_current]=Rectangular_QAM_slicer(symOut_current, M, pav);
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
for i=1:Nt
    symOut_tmp=symOut_prev;
    symOut_tmp(i)=0;
    y_tmp=y-H*W*symOut_tmp;
    symOut_current(i)=((H(:,i)')/(H(:,i)'*H(:,i)+SNRd^(-1)))*y_tmp;
%     symI(i)=sym2(i);
end
[symOut_current]=Rectangular_QAM_slicer(symOut_current,M, pav);
%   symOut=sym2;
%   return;
if(stage==maxStage)
    symOut=symOut_current;
    return
end
[symOut]=MIC_Recursive(y, H, symOut_current, SNRd, M, pav, stage+1, maxStage, W,G, symConstell);
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

function [X_hat] = Rectangular_QAM_slicer(X, M, pav)    %slicer
%need to be modified

% sq10=sqrt(10);
if(M==2)
    d=sqrt(pav);
else
d=sqrt(3*pav/(2*(M-1)));
end

N=length(X);
for i=1:N
if(M==2)
    
    if(real(X(i))<=0)
        REAL=-d;
    else
        REAL=d;
    end
    X(i)=complex(REAL,0);
elseif(M==4)
        if(real(X(i))<=0)
        REAL=-d;
    else
        REAL=d;
        end
        if(imag(X(i))<=0)
        IMAG=-d;
    else
        IMAG=d;
        end
        X(i)=complex(REAL, IMAG);
elseif(M==16)
    
        if(real(X(i))<=-2*d)
        REAL=-3*d;
        elseif (real(X(i))>-2*d&&real(X(i))<=0)
        REAL=-d;
        elseif(real(X(i))>0&&real(X(i))<=2*d)
            REAL=d;
        elseif(real(X(i)>2*d))
            REAL=3*d;
        end
        if(imag(X(i))<=-2*d)
        IMAG=-3*d;
        elseif (imag(X(i))>-2*d&&imag(X(i))<=0)
        IMAG=-d;
        elseif(imag(X(i))>0&&imag(X(i))<=2*d)
            IMAG=d;
        elseif(imag(X(i)>2*d))
            IMAG=3*d;
        end
        X(i)=complex(REAL, IMAG);
end

end
X_hat=X;
end



