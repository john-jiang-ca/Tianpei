function [ symOut] = RSVR( H,  y, SNRd,  M, pav)
%real SVR detection 
%   Output
% symOut: regression coefficiences estimation
% lamida: dual variable
% Theta: the vector that store all the objective function value
% G: duality Gap
% iteration: iteration time
% MSE: mean square error
%Input
% pH: data sample
% y: data  sample output
% SNRd: signal to noise ratio at output 
% symbolConstellation: the symbol constellation
% iteration: iteration time
% epsilon=1e-7;   %epsilon
% tol=1e-2;  %tolerance of duality gap
tol=1e-2;
C=4;
epsilon=1e-7;
K=H*H';
Nt=length(H(1,:));
Nr=length(K(1,:));
% C=5;
Phi=zeros(Nr,1);
sigma=zeros(Nr,1);
lambda=zeros(Nr,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
Theta(1)=0;
Phi_tmp=abs(Phi)-epsilon;
Phi_tmp=Phi_tmp(find(Phi_tmp>0));
Chi=norm(Phi_tmp)^(2);
Theta(1)=-0.5*C*Chi;
G(1)=-2*Theta(1);

% G(1)=G(1)/abs(G(1)+Theta(1));
%% update lamida

iteration=1;
while(G(iteration)>tol)

delta=zeros(Nr,1);
sigma_tmp=zeros(Nr, 1);
NoiseTerm=0;
sgn=zeros(Nr,1);
K_tmp=0;
 %% approximate 2-D work set selection (without damping-bubble sorting)
for count=1:Nr
%      k=(K(count,count)+1/C);
%        lamida_tmp1=lamida(count)*(K(count,count)/k)+(Phi(count)-epsilon*(-1))/k;
%       lamida_tmp2=lamida(count)*(K(count,count)/k)+(Phi(count)-epsilon*(1))/k;
  lambda_tmp1=lambda(count)+(Phi(count)-epsilon*(-1))/(K(count,count));
  lamida_tmp2=lambda(count)+(Phi(count)-epsilon*(1))/(K(count,count));
        %clipping
%     if(lamida_tmp1<-C)
%         lamida_tmp1=-C;
%     elseif(lamida_tmp1>C)
%         lamida_tmp1=C;
%     end
        %clipping
%     if(lamida_tmp2<-C)
%         lamida_tmp2=-C;
%     elseif(lamida_tmp2>C)
%         lamida_tmp2=C;
%     end
    delta1=(lambda_tmp1-lambda(count))*((-1/2)*(lambda_tmp1-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lambda_tmp1)-abs(lambda(count)));
%     -(1/C)*(lamida_tmp1^2-(lamida(count)^2));
    delta2=(lamida_tmp2-lambda(count))*((-1/2)*(lamida_tmp2-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp2)-abs(lambda(count)));
%     -(1/C)*(lamida_tmp1^2-(lamida(count)^2));
    if(delta1>delta2)
        delta(count)=delta1;
        sgn(count)=-1;
%         lamida_tmp=lamida_tmp1;
        
    else
        delta(count)=delta2;
        sgn(count)=1;
%         lamida_tmp=lamida_tmp2;
    end

%     sigma_tmp(count)=lamida_tmp-lamida(count);
end
%bubble sorting
%find the coordinates that can generate first and second largest
%objective gain
[value, First]=max(delta);
delta(First)=min(delta)-1;
[value, Second]=max(delta);


    


%% update lamida, sigma, Phi, Theta and duality gap
K_tmp=K(First,First)*K(Second,Second)-(K(First,Second))^2;
%update lamida(First)
lambda_tmp1=lambda(First)+(Phi(First)*K(Second,Second)-Phi(Second)...
    *K(First,Second)-epsilon*(sgn(First)*K(Second,Second)...
    -sgn(Second)*K(First,Second)))/K_tmp;
sigma(First)=lambda_tmp1-lambda(First);
lambda(First)=lambda_tmp1;

%update lamida(Second)
lamida_tmp2=lambda(Second)+(Phi(Second)*K(First,First)-Phi(First)...
    *K(First,Second)-epsilon*(sgn(Second)*K(First,First)...
    -sgn(First)*K(First,Second)))/K_tmp;
sigma(Second)=lamida_tmp2-lambda(Second);
lambda(Second)=lamida_tmp2;
% for count=1:Nr
%     Phi(count)=Phi(count)-sigma(First)*K(count,First)-sigma(Second)*K(count,Second);
% end
Phi=Phi-sigma(First)*K(:,First)-sigma(Second)*K(:,Second);
iteration=iteration+1;
% for count1=1:Nr
%      if(abs(Phi(count1))-epsilon>0)
%      NoiseTerm=NoiseTerm+(abs(Phi(count1))-epsilon)^2;
%      end
% % NoiseTerm=NoiseTerm+(lamida(count))^2;
% end
Chi_tmp=Phi-epsilon;
Chi_tmp=Chi_tmp(find(Chi_tmp>0));
NoiseTerm=norm(Chi_tmp)^(2);
% NoiseTerm=1/C*NoiseTerm;
%calculate objective function value
Theta_tmp=-(1/2)*lambda'*K*lambda+y'*lambda-epsilon*norm(lambda, 1);
Theta(iteration)=Theta_tmp-0.5*C*NoiseTerm;
%calculate duality gap
G(iteration)=y'*lambda-epsilon*(norm(lambda,1))-2*Theta(iteration);
G(iteration)=G(iteration)/(G(iteration)+Theta(iteration));



end
symOut_tmp=H'*lambda;
if(length(symOut_tmp)==2)
    symOut=symOut_tmp(1)+1i*symOut_tmp(2);
else
symOut=symOut_tmp(1:Nt/2)+1i*symOut_tmp(Nt/2+1:Nt);
end
symOut=Rectangular_QAM_slicer(symOut, M, pav);
% MSE=norm(y-H*symOut);
%   symOut_MMSE=zeros(Nt,1);
%   I=ones(Nt,Nt);
%   symOut_MMSE=(inv(pH'*pH)+SNRd^(-1)*I)*pH'*y;
%% rounding 


%% mean square error
% MSE=norm(y-pH*symOut);
end
