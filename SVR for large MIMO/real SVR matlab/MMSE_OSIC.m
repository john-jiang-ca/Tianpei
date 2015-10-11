function [ symOut ] = MMSE_OSIC(y, H, SNRd, M, pav)
%this routine implement MMSE-OSIC algorithm
%   Detailed explanation goes here

Nr=length(H);
Nt=length(H(:,1));
symOut=zeros(Nt,1);  %output symbol vector
% symOut_tmp_v=zeros(Nt,1);  %
symOut_list=1:Nt; %record the current index of the updated subvector
% symOut_Recons=zeros(Nt,1);

H_tmp=H;  %updated channel matrix
y_tmp=y; %updated received symbol vector 
for count=1:Nt
    I=eye(Nt-(count-1));
    G=(I/(H_tmp'*H_tmp+SNRd^(-1)*I));   %inverse of HHH matrix with SNR
    [value, index]=min(diag(G));   %choice of the channel with the minimum post processing SNR
    G=G*H_tmp';   %equalization matrix
    symOut_tmp=G(index,:)*y_tmp;   
    symOut_tmp=Rectangular_QAM_slicer(symOut_tmp,M,pav); %detected symbol 
    symOut(symOut_list(index))=symOut_tmp; %upate output symbol vector
    y_tmp=y_tmp-H_tmp(:,index)*symOut_tmp; %update received symbol vector
    H_tmp(:,index)=[];  %update channel matrix
    symOut_list(index)=[];   %update symbol index list

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

