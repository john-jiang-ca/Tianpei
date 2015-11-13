function [ symOut ] = CSVD( y, H, SNRd, M, pav, C, epsilon, tol )
%omplex support vector detector
%   Detailed explanation goes here
%% hyperparameter settings
% epsilon=1e-2;
% tol=1e-2;
% C=1;
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
Theta(1)=(-0.5)*C*(ChiR+ChiI);
G(1)=-2*Theta;
ratio(1)=abs(G(1))/(abs(G(1))+abs(Theta(1)));
iteration=1;
while (ratio(iteration)>=tol&& iteration<=1e3)
%% double channel solver
% LambdaR_new=zeros(Nr, 1);
% LambdaI_new=zeros(Nr,1);
% PhiR_new=zeros(Nr,1);
% PhiI_new=zeros(Nr,1);
% PsiR_new=zeros(Nr,1);
% PsiI_new=zeros(Nr,1);
% i=0;
% j=0;
% m=0;
% n=0;
%% real channel
labelR=1;
sigmaR1=0;
sigmaR2=0;
sigmaR1_a=0;
sigmaR2_a=0;
[i,j, sigmaR1, sigmaR2, sigmaR1_a, sigmaR2_a, LambdaR_new, PhiR_new, PsiI_new, DeltaThetaR]...
    =a2DSolver(Kr, Ki, LambdaR, PhiR, PsiI, epsilon, labelR);
% [index1,index2, sigma1, sigma2, sigma1_a, sigma2_a, Lambda_new, Phi_new, Psi_new, deltaTheta]=a2DSolver(Kr, Ki, Lambda, Phi, Psi, epsilon,  label)
%% imaginary channel 
labelI=0;
sigmaI1=0;
sigmI2=0;
sigmaI1_a=0;
sigmaI2_a=0;
[m,n, sigmaI1, sigmaI2, sigmaI1_a, sigmaI2_a, LambdaI_new, PhiI_new, PsiR_new, DeltaThetaI]...
    =a2DSolver(Kr, Ki, LambdaI, PhiI, PsiR, epsilon, labelI);





%% update G Theta and Phi Psi
Theta_testR(iteration)=DeltaThetaR;
Theta_testI(iteration)=DeltaThetaI;

%% real channel update
ChiR_tmp=abs(PhiR_new+PsiR_new)-epsilon;   %update noise term
ChiR_tmp=ChiR_tmp(find(ChiR_tmp>0));
ChiR_new=norm(ChiR_tmp)^(2);
ChiR_v(iteration)=ChiR_new;
Delta_ChiR=ChiR_new-ChiR;    %noise term change of real channel 
Delta_ThetaR=DeltaThetaR-0.5*C*Delta_ChiR;    %Theta change of real channel
DeltaThetaR_v(iteration)=Delta_ThetaR;
Delta_GR=real(y(i))*sigmaR1+real(y(j))*sigmaR2-epsilon*(sigmaR1_a+sigmaR2_a)-2*Delta_ThetaR;   %change  of G in real channel 
Delta_GR_v(iteration)=Delta_GR;


%% imaginary channel update
ChiI_tmp=abs(PhiI_new+PsiI_new)-epsilon;
ChiI_tmp=ChiI_tmp(find(ChiI_tmp>0));
ChiI_new=norm(ChiI_tmp)^(2);
ChiI_v(iteration)=ChiI_new;
Delta_ChiI=ChiI_new-ChiI;    %noise term change of imaginary channel
Delta_ThetaI=DeltaThetaI-0.5*C*Delta_ChiI;   %Theta change of imaginary channel
DeltaThetaI_v(iteration)=Delta_ThetaI;
Delta_GI=imag(y(m))*sigmaI1+imag(y(n))*sigmaI2-epsilon*(sigmaI1_a+sigmaI2_a)-2*Delta_ThetaI;   %change of G in imaginary channel
Delta_GI_v(iteration)=Delta_GI;

%% total update
Delta_Chi=Delta_ChiR+Delta_ChiI;   %total  noise term change
delta_Chi_v(iteration)=Delta_Chi;
Delta_Theta=Delta_ThetaR+Delta_ThetaI; %calculate the update value of the total Theta
Delta_G=Delta_GR+Delta_GI;    %total G change
DeltaG_v(iteration)=Delta_G;

%calculate real part of K correlated G
DeltaG_add=sigmaR1*Ki(i,m)*sigmaI1+sigmaR2*Ki(j,m)*sigmaI1+sigmaR1*Ki(i,n)*sigmaI2...
    +sigmaR2*Ki(j,n)*sigmaI2+sigmaR1*PsiR(i)+sigmaR2*PsiR(j)+sigmaI1*PsiI(m)+sigmaI2*PsiI(n);   %addtional change in total G





%% update final G and Theta
Delta_G=Delta_G-DeltaG_add;
iteration=iteration+1;
Theta(iteration)=Theta(iteration-1)+Delta_Theta;
G(iteration)=G(iteration-1)+Delta_G;
% G_test(iteration)=G(iteration)+DeltaG_add;
%calculate imaginary part of K correlated G
% deltaG=deltaG+deltaG_additional;
% Theta=Theta+DeltaTheta;
% G=G+deltaG;
% ratio=G/(G+Theta);
%% update Phi and Psi of real and imaginary channel 
PhiR=PhiR_new;
PsiR=PsiR_new;
PhiI=PhiI_new;
PsiI=PsiI_new;
LambdaR=LambdaR_new;
LambdaI=LambdaI_new;
% Theta_tmp=(-1/2)*LambdaR'*Kr*LambdaR+(-1/2)*LambdaI'*Kr*LambdaI+real(y)'*LambdaR+imag(y)'*LambdaI...
%     -epsilon*(norm(LambdaR,1)+norm(LambdaI,1))-(C/2)*(ChiR_new+ChiI_new); %check the correctness of the updata rule for G and Theta
ChiR=ChiR_new;
ChiI=ChiI_new;

% G_tmp=real(y)'*LambdaR+imag(y)'*LambdaI-epsilon*(norm(LambdaR,1)+norm(LambdaI,1))-2*Theta_tmp-LambdaR'*Ki*LambdaI; 

%% calculate ratio
if(ratio(iteration-1)<1e-2)
    controller=1e-2;
else
    controller=1;
end
ratio(iteration)=controller*abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));
% iteration=iteration+1;
% G_v=[G_v,Theta1];
end
%% reconstruct output symbol vector
Lambda=complex(LambdaR, LambdaI);
x_hat=H'*(Lambda);
[symOut]=Rectangular_QAM_slicer(x_hat, M, pav);
end



