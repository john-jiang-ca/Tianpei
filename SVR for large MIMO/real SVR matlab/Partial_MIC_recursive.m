function [ symOut ] = Partial_MIC_recursive( y, H, symOut_total, list, H1, symOut_prev, SNRd, M, pav, stage, maxStage, symConstell, tol)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%Multistage Interference Cancellation (MIC) recursive algorithm
%   Detailed explanation goes here

Nt=length(H(1,:)); %number of receive antennas
% Nt=length(H(:,1));  %number of transmit antennas
I=eye(Nt);
%% IC
if(stage~=1)
y=y-H1*symOut_prev;
end
%% MMSE
symOut_tmp=(I/(H'*H+SNRd^(-1)*I))*H'*y;
[symOut_tmp1]=Rectangular_QAM_slicer(symOut_tmp,M, pav);
%% check reliability
[W]=weightCal(symOut_tmp, symOut_tmp1, tol, pav, M);
%% update H and symOut_total
count1=1;
symOut_un=zeros(length(symOut_total),1);
H_tmp=[];
list_tmp=[];
symOut_current=[];
H1=[];
for count=1:Nt
    if(W(count)==1)
        H1=[H1,H(:,count)];
%         H(:,count)=[];
        symOut_total(list(count))=symOut_tmp1(count);
        symOut_current=[symOut_current; symOut_tmp1(count)];
%         list(count)=[];
        count1=count1+1;
    else
    symOut_un(list(count))=symOut_tmp1(count);
    H_tmp=[H_tmp,H(:,count)];
    list_tmp=[list_tmp, list(count)];
    end
end

% maxStage=4;
%% First stage estimation by MMSE
% if(Nstage==0)
% symOut_current=zeros(Nt,1);
% if(stage==1)
% INV=(I/(H'*H+SNRd^(-1)*I));
% G=INV*H';
% symOut_current=G*y;
% [symOut_current]=Rectangular_QAM_slicer(symOut_current, M, pav);
% symOut_prev=symOut_current;
% end


if(stage==maxStage)
    if(~isempty(list))
        for count=1:length(list)
    symOut_total(list(count))=symOut_un(list(count));
        end
    end
    symOut=symOut_total;
    return
end
[symOut]=Partial_MIC_recursive(y, H_tmp, symOut_total, list_tmp, H1, symOut_current, SNRd, M, pav, stage+1, maxStage, symConstell, tol);
end


%% adaptive  weight function 
function [W]=weightCal(s, s_Q, tol, pav, M)
 d=sqrt(3*pav/(2*(M-1)));
 L=length(s);
 W=zeros(L,1);
 for count=1:L
     if(abs(real(s(count))-real(s_Q(count)))>tol*d||abs(imag(s(count))-imag(s_Q(count)))>tol*d)
         W(count)=0;
     else
         W(count)=1;
     end
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


