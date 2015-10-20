function [ symOut ] = BI_GDFE( y, H, SNRd, M, pav )
%this routine implement block iterative generalized decision feedback
%equalizer
%   Detailed explanation goes here
Nt=length(H(1,:));
Nr=length(H(:,1));
maxIteration=4;   % the maximal iteration time
p=[0.95,0.95,0.99,0.999];  % the correlation factor of the estimated symbol vector of ith iteration with the true symbol vector
symOut_m=zeros(Nt, maxIteration);   % the matrix to store all the estimations in all iterations
%% get preliminary estimation
I=eye(Nr);
G=(I/(H'*H+SNRd^(-1)*I))*H';
symOut_pre=G*y;
symOut_pre=Rectangular_QAM_slicer(symOut_pre, M, pav);

%% block iteration with FFE and FBE
symOut_prev=symOut_pre;
for count=1:maxIteration
    K=(I/((1-p(count)^(2))*H'*H+SNRd^(-1)*I))*H';
    D=p(count)*(diag(K*H)-K*H);
    symOut_m(:,count)=K*y+D*symOut_prev;
    symOut_m(:,count)=Rectangular_QAM_slicer(symOut_m(:,count), M, pav);
    symOut_prev=symOut_m(:,count);
end

%% the final decision 
symOut=symOut_prev;
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