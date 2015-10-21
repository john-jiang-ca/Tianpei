function [ symOut ] = CSVD( y, H, SNRd, M, G_v )
%omplex support vector detector
%   Detailed explanation goes here
%% hyperparameter settings
epsilon=1e-2;
tol=1e-2;
C=1;
G_v=[];
%% initialization 
Nr=length(H(1,:));
Nt=length(H(:,1));
K=H*H';
Kr=real(K);
Ki=imag(K);
Lambda=zeros(Nr,1);
LambdaR=zeros(Nr,1);
LambdaI=zeros(Nr,1);
Phi=zeros(Nr,1);
PhiR=real(y);
PhiI=imag(y);
Phi=PhiR+1i*PhiI;
Psi=zeros(Nr,1);
PsiR=zeros(Nr,1);
PsiI=zeros(Nr,1);
Psi=PhiR+1i*PsiI;
% ChiR=0;
% ChiI=0;
ChiR_tmp=abs(PhiR+PsiR)-epsilon;
ChiR_tmp=ChiR_tmp(find(ChiR_tmp>0));
ChiR=norm(ChiR_tmp, 2)^(2);
ChiI_tmp=abs(PhiI+PsiI)-epsilon;
ChiI_tmp=ChiI_tmp(find(ChiI_tmp>0));
ChiI=norm(ChiI_tmp,2)^(2);
Theta=(-C/2)*(ChiR+ChiI);
G=-2*Theta;
ratio=G/(G+Theta);
iteration=0;
while (ratio>=tol&&iteration<1e3)
%% double channel solver
LambdaR_new=zeros(Nr, 1);
LambdaI_new=zeros(Nr,1);
PhiR_new=zeros(Nr,1);
PhiI_new=zeros(Nr,1);
PsiR_new=zeros(Nr,1);
PsiI_new=zeros(Nr,1);
i=0;
j=0;
m=0;
n=0;
%% real channel
labelR=1;
sigmaR1=0;
sigmaR2=0;
sigmaR1_a=0;
sigmaR2_a=0;
[i,j, sigmaR1, sigmaR2, sigmaR1_a, sigmaR2_a, LambdaR_new, PhiR_new, PsiI_new, DeltaThetaR]=a2DSolver(Kr, Ki, LambdaR, PhiR, PsiI, epsilon, labelR);
% [index1,index2, sigma1, sigma2, sigma1_a, sigma2_a, Lambda_new, Phi_new, Psi_new, deltaTheta]=a2DSolver(Kr, Ki, Lambda, Phi, Psi, epsilon,  label)
%% imaginary channel 
labelI=0;
sigmaI1=0;
sigmI2=0;
sigmaI1_a=0;
sigmaI2_a=0;
[m,n, sigmaI1, sigmaI2, sigmaI1_a, sigmaI2_a, LambdaI_new, PhiI_new, PsiR_new, DeltaThetaI]=a2DSolver(Kr, Ki, LambdaI, PhiI, PsiR, epsilon, labelI);
%% update G and Theta
% sigmaR1=LambdaR_new(i)-LambdaR(i);  %calculate the sigma=lambda^{new}_{i}-lambda_{i}
% sigmaR1_a=abs(LambdaR_new(i))-abs(LambdaR(i));
% sigmaR2=LambdaR_new(j)-LambdaR(j);
% sigmaR2_a=abs(LambdaR_new(j))-abs(LambdaR(j));
% sigmaI1=LambdaI_new(m)-LambdaI(m);
% sigmaI1_a=abs(LambdaI_new(m))-abs(LambdaI(m));
% sigmaI2=LambdaI_new(n)-LambdaI(n);
% sigmaI2_a=abs(LambdaI_new(n))-abs(LambdaI(n));
ChiR_tmp=abs(PhiR_new+PsiR_new)-epsilon;   %update noise term
ChiR_tmp=ChiR_tmp(find(ChiR_tmp>0));
ChiR_new=norm(ChiR_tmp,2)^(2);
ChiI_tmp=abs(PhiI_new+PsiI_new)-epsilon;
ChiI_tmp=ChiI_tmp(find(ChiI_tmp>0));
ChiI_new=norm(ChiI_tmp,2)^(2);
DeltaTheta=DeltaThetaR+DeltaThetaI-(C/2)*(ChiR_new+ChiI_new-ChiR-ChiI); %calculate the update value of delta Theta
deltaG=real(y(i))*sigmaR1+real(y(j))*sigmaR2+imag(y(m))*sigmaI1+imag(y(n))*sigmaI2-epsilon*(sigmaR1_a+sigmaR2_a+sigmaI1_a+sigmaI2_a)-2*DeltaTheta; 
%calculate real part of K correlated G
deltaG_additional=sigmaR1*Ki(i,m)*sigmaI1+sigmaR2*Ki(j,m)*sigmaI1+sigmaR1*Ki(i,n)*sigmaI2+sigmaR2*Ki(j,n)*sigmaI2+sigmaR1*PsiR(i)+sigmaR2*PsiR(j)+sigmaI1*PsiI(m)+sigmaI2*PsiI(n);
%calculate imaginary part of K correlated G
% deltaG=deltaG+deltaG_additional;
% Theta=Theta+DeltaTheta;
% G=G+deltaG;
% ratio=G/(G+Theta);
PhiR=PhiR_new;
PsiR=PsiR_new;
PhiI=PhiI_new;
PsiI=PsiI_new;
LambdaR=LambdaR_new;
LambdaI=LambdaI_new;
ChiR=ChiR_new;
ChiI=ChiI_new;
Theta1=(-1/2)*LambdaR'*Kr*LambdaR+(-1/2)*LambdaI'*Kr*LambdaI+real(y)'*LambdaR+imag(y)'*LambdaI-epsilon*(norm(LambdaR,1)+norm(LambdaI,1))-(C/2)*(ChiR+ChiI);
G1=real(y)'*LambdaR+imag(y)'*LambdaI-epsilon*(norm(LambdaR,1)+norm(LambdaI,1))-2*Theta-LambdaR'*Ki*LambdaI;
ratio=G1/(G1+Theta1);
iteration=iteration+1;
G_v=[G_v,Theta1];
end
%% reconstruct output symbol vector
Lambda=LambdaR+1i*LambdaI;
x_hat=H'*(Lambda);
[symOut]=Rectangular_QAM_slicer(x_hat, M, 1/Nt);
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

