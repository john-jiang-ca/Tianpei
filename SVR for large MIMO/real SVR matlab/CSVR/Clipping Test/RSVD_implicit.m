function [ symOut] = RSVD_implicit( H,  y, SNRd,  M, pav, C, tol , epsilon)
%real SVR detection with uniform formula 
%noise term is considered uniformly in the dual objective function as in Smola's paper 
%   Output
% symOut: regression coefficiences estimation
% lambda: dual variable
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

K=H*H';

N=length(H(1,:));
L=length(K(1,:));
I=eye(L);
K_hat=K+(1/C)*I;
% C=5;
Phi=zeros(L,1);
sigma=zeros(L,1);
lambda=zeros(L,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
Theta(1)=0;
% Phi_tmp=abs(Phi)-epsilon;
% Phi_tmp=Phi_tmp(find(Phi_tmp>0));
% Chi(1)=norm(Phi_tmp)^(2);
% Theta(1)=-0.5*C*Chi(1);
% G(1)=-2*Theta(1);

% G(1)=G(1)/abs(G(1)+Theta(1));
%% update lambda
% ratio(1)=G(1)/(G(1)+Theta(1));
ratio(1)=1;
iteration=1;
while(ratio(iteration)>tol)

delta=zeros(L,1);
sigma_tmp=zeros(L, 1);
NoiseTerm=0;
sgn=zeros(L,1);
K_tmp=0;
 %% approximate 2-D work set selection (without damping-bubble sorting)
for count=1:L
%      k=(K(count,count)+1/C);
%        lambda_tmp1=lambda(count)*(K(count,count)/k)+(Phi(count)-epsilon*(-1))/k;
%       lambda_tmp2=lambda(count)*(K(count,count)/k)+(Phi(count)-epsilon*(1))/k;
  lambda_tmp1=lambda(count)*K(count,count)/K_hat(count,count)+(Phi(count)-epsilon*(-1))/(K_hat(count,count));
  lambda_tmp2=lambda(count)*K(count,count)/K_hat(count,count)+(Phi(count)-epsilon*(1))/(K_hat(count,count));
% %         clipping
% xi_tmp_upper=max(0, Phi(count)-epsilon);
% xi_tmp_lower=max(0, -Phi(count)-epsilon);
% clipper_upper=C*xi_tmp_upper;
% clipper_lower=-C*xi_tmp_lower;
%     if(lambda_tmp1<-clipper_lower)
%         lambda_tmp1=-clipper_lower;
%     elseif(lambda_tmp1>clipper_upper)
%         lambda_tmp1=clipper_upper;
%     end
% % clipping
%     if(lambda_tmp2<-clipper_lower)
%         lambda_tmp2=-clipper_lower;
%     elseif(lambda_tmp2>clipper_upper)
%         lambda_tmp2=clipper_upper;
%     end
    delta1=(lambda_tmp1-lambda(count))*((-1/2)*(lambda_tmp1-lambda(count))*K_hat(count,count)+Phi(count)-(1/C)*lambda_tmp1)-epsilon*(abs(lambda_tmp1)...
        -abs(lambda(count)));
%     -(1/C)*(lambda_tmp1^2-(lambda(count)^2));
    delta2=(lambda_tmp2-lambda(count))*((-1/2)*(lambda_tmp2-lambda(count))*K_hat(count,count)+Phi(count)-(1/C)*lambda_tmp2)-epsilon*(abs(lambda_tmp2)...
        -abs(lambda(count)));
%     -(1/C)*(lambda_tmp1^2-(lambda(count)^2));
    if(delta1>delta2)
        delta(count)=delta1;
        sgn(count)=-1;
%         lambda_tmp=lambda_tmp1;
        
    else
        delta(count)=delta2;
        sgn(count)=1;
%         lambda_tmp=lambda_tmp2;
    end

%     sigma_tmp(count)=lambda_tmp-lambda(count);
end
%bubble sorting
%find the coordinates that can generate first and second largest
%objective gain
[value, First]=max(delta);
delta(First)=min(delta)-1;
[value, Second]=max(delta);


    


%% update lambda, sigma, Phi, Theta and duality gap
K_tmp=K_hat(First,First)*K_hat(Second,Second)-(K(First,Second))^2;
M1=Phi(First)-epsilon*sgn(First);
M2=Phi(Second)-epsilon*sgn(Second);
%update lambda(First)
% lambda_tmp1=lambda(First)+(Phi(First)*K(Second,Second)-Phi(Second)...
%     *K(First,Second)-epsilon*(sgn(First)*K(Second,Second)-sgn(Second)*K(First,Second)))/K_tmp;
lambda_tmp1=(K(First, First)*K_hat(Second, Second)-K(First, Second)^(2))*lambda(First)+lambda(Second)*K(First, Second)...
    *(K_hat(Second, Second)-K(Second, Second))+M1*K_hat(Second,Second)-M2*K(First, Second);
lambda_tmp1=lambda_tmp1/K_tmp;
sigma(First)=lambda_tmp1-lambda(First);
sigma1_abs=abs(lambda_tmp1)-abs(lambda(First));
lambda(First)=lambda_tmp1;

%update lambda(Second)
% lambda_tmp2=lambda(Second)+(Phi(Second)*K(First,First)-Phi(First)...
%     *K(First,Second)-epsilon*(sgn(Second)*K(First,First)...
%     -sgn(First)*K(First,Second)))/K_tmp;
lambda_tmp2=(K(Second, Second)*K_hat(First, First)-K(First, Second)^(2))*lambda(Second)+lambda(First)*K(First, Second)...
    *(K_hat(First, First)-K(First, First))+M2*K_hat(First,First)-M1*K(First, Second);
lambda_tmp2=lambda_tmp2/K_tmp;
sigma(Second)=lambda_tmp2-lambda(Second);
sigma2_abs=abs(lambda_tmp2)-abs(lambda(Second));
lambda(Second)=lambda_tmp2;
% for count=1:Nr
%     Phi(count)=Phi(count)-sigma(First)*K(count,First)-sigma(Second)*K(count,Second);
% end

iteration=iteration+1;
% for count1=1:Nr
%      if(abs(Phi(count1))-epsilon>0)
%      NoiseTerm=NoiseTerm+(abs(Phi(count1))-epsilon)^2;
%      end
% % NoiseTerm=NoiseTerm+(lambda(count))^2;
% end
%% update Theta G and Phi
Phi_new=Phi-sigma(First)*K(:,First)-sigma(Second)*K(:,Second);
% Chi_tmp=abs(Phi_new)-epsilon;
% Chi_tmp=Chi_tmp(find(Chi_tmp>0));
% Chi(iteration)=norm(Chi_tmp)^(2);
% Delta_Chi=Chi(iteration)-Chi(iteration-1);
% NoiseTerm=1/C*NoiseTerm;
%calculate objective function value
% Theta_tmp1=-(1/2)*lambda'*K*lambda+y'*lambda-epsilon*norm(lambda, 1)-0.5*C*Chi(iteration);
% Delta_Theta=(-1/2)*(sigma(First)^(2)*K_hat(First, First)+sigma(Second)^(2)*K_hat(Second,Second)+2*sigma(First)*sigma(Second)*K(First,Second))...
%     +(Phi(First)-(1/C)*Lambda(First))*sigma(First)+sigma(Second)*(Phi(Second)-(1/C)*Lambda(Second))-epsilon*(sigma1_abs+sigma2_abs);
% Delta_Theta=-(1/2)*(K_hat(First, First)*sigma(First)^(2)+K_hat(Second, Second)*sigma(Second)^(2)+2*K(First,Second)*sigma(First)*sigma(Second))...
%     +(Phi(First)-(1/C)*lambda(First))*sigma(First)+(Phi(Second)-(1/C)*lambda(Second))*sigma(Second)-epsilon*(sigma1_abs+sigma2_abs);
% Delta_Theta=Delta_Theta-0.5*C*Delta_Chi;
% Theta_tmp=-(1/2)*lambda'*K_hat*lambda+y'*lambda-epsilon*norm(lambda,1);
% Theta(iteration)=Theta(iteration-1)+Delta_Theta;
% %calculate duality gap
% %  G_tmp=y'*lambda-epsilon*(norm(lambda,1))-2*Theta(iteration);
% Delta_G=y(First)*sigma(First)+y(Second)*sigma(Second)-epsilon*(sigma1_abs+sigma2_abs)-2*Delta_Theta;
% G_tmp=lambda'*K_hat*lambda-y'*lambda+epsilon*(norm(lambda,1))+(1/C)*(lambda'*lambda);
% G(iteration)=G(iteration-1)+Delta_G;

Phi_new=Phi-sigma(First)*K(:,First)-sigma(Second)*K(:,Second);
Chi_tmp=abs(Phi_new)-epsilon;
Chi_tmp=Chi_tmp(find(Chi_tmp>0));
Chi(iteration)=norm(Chi_tmp)^(2);
Delta_Chi=Chi(iteration)-Chi(iteration-1);  %the updated value of Chi
Delta_Theta=(-1/2)*(sigma(First)^(2)*K(First, First)+sigma(Second)^(2)*K(Second,Second)+2*sigma(First)*sigma(Second)*K(First,Second))...
    +Phi(First)*sigma(First)+sigma(Second)*Phi(Second)-epsilon*(sigma1_abs+sigma2_abs);  
Delta_Theta=Delta_Theta-0.5*C*Delta_Chi;   %the gain of Theta
Theta(iteration)=Theta(iteration-1)+Delta_Theta;  %update Theta
Delta_G=y(First)*sigma(First)+y(Second)*sigma(Second)-epsilon*(sigma1_abs+sigma2_abs)-2*Delta_Theta; %gain of duality gap
G(iteration)=G(iteration-1)+Delta_G;  %update G
if(ratio(iteration-1)<1e-2)
controller=1e-2;
else
    controller=1;
end
% ratio(iteration)=controller*abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));
ratio(iteration)=abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));
Phi=Phi_new;


end
%% Reconstruct regression coefficience vector based on lambda
symOut_tmp=H'*lambda;
if(length(symOut_tmp)==2)
    symOut=symOut_tmp(1)+1i*symOut_tmp(2);
else
symOut=symOut_tmp(1:N/2)+1i*symOut_tmp(N/2+1:N);
end
symOut=Rectangular_QAM_slicer(symOut, M, pav);   %slicing
end