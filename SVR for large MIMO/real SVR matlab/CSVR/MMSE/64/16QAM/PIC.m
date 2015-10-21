function [ symOut ] = PIC( y, H, SNRd,M, pav)
%This function performs original parallel interference cancellation (PIC)
%algorithm
%   Detailed explanation goes here
% Nstage=0;  %number of stages
Nr=length(H); %number of receive antennas
Nt=length(H(:,1));  %number of transmit antennas
I=eye(Nt);
%% First stage estimation by MMSE
% if(Nstage==0)
symI=(I/(H'*H+SNRd^(-1)*I))*H'*y;
[symI]=Rectangular_QAM_slicer(symI,M, pav);
%    symOut=symI;
%    return;
% Nstage=Nstage+1;
% end
% if(Nstage==MAXstage)
%     return;
% end
%% 2-stage PIC
sym2=zeros(Nt,1);
for i=1:Nt
    symI_tmp=symI;
    symI_tmp(i)=0;
    y_tmp=y-H*symI_tmp;
    sym2(i)=(H(:,i)'/(H(:,i)'*H(:,i)+SNRd^(-1)))*y_tmp;
%     symI(i)=sym2(i);
end
[sym2]=Rectangular_QAM_slicer(sym2,M, pav);
%   symOut=sym2;
%   return;
%% 3-stage PIC
sym3=zeros(Nt,1);
for i=1:Nt
    symI_tmp=sym2;
    symI_tmp(i)=0;
    y_tmp=y-H*symI_tmp;
    sym3(i)=(H(:,i)'/(H(:,i)'*H(:,i)+SNRd^(-1)))*y_tmp;
%     sym2(i)=sym3(i);
end
 [sym3]=Rectangular_QAM_slicer(sym3,M, pav);
% [symI]=PIC(y,H,SNRd, Nstage+1,MAXstage);
%% 4-stage PIC 
sym4=zeros(Nt,1);
for i=1:Nt
    symI_tmp=sym3;
    symI_tmp(i)=0;
    y_tmp=y-H*symI_tmp;
    sym4(i)=(H(:,i)'/(H(:,i)'*H(:,i)+SNRd^(-1)))*y_tmp;
%     sym3(i)=sym3(i);
end
[symOut]=Rectangular_QAM_slicer(sym3,M, pav);
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
