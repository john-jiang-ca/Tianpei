function [ W_e ] = SVR_learning(y ,x , SNRd,C ,epsilon, tol )
%This subroutine implement the Support Vector regression 
%   Detailed explanation goes here
y=[real(y); imag(y)];
x=[real(x), -imag(x); imag(x), real(x)];
N=length(x(1,:));
K=x*x';
% Nt=length(x(1,:));
M=length(K(1,:));
% C=5;
Phi=zeros(M,1);
sigma=zeros(M,1);
lambda=zeros(M,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
Theta(1)=0;
Phi_tmp=abs(Phi)-epsilon;
Phi_tmp=Phi_tmp(find(Phi_tmp>0));
Chi(1)=norm(Phi_tmp)^(2);
Theta(1)=-0.5*C*Chi(1);
G(1)=-2*Theta(1);

% G(1)=G(1)/abs(G(1)+Theta(1));
%% update lambda
ratio(1)=G(1)/(G(1)+Theta(1));
iteration=1;
while(ratio(iteration)>tol)

delta=zeros(M,1);
sigma_tmp=zeros(M, 1);
NoiseTerm=0;
sgn=zeros(M,1);
K_tmp=0;
 %% approximate 2-D work set selection (without damping-bubble sorting)
for count=1:M
%      k=(K(count,count)+1/C);
%        lambda_tmp1=lambda(count)*(K(count,count)/k)+(Phi(count)-epsilon*(-1))/k;
%       lambda_tmp2=lambda(count)*(K(count,count)/k)+(Phi(count)-epsilon*(1))/k;
  lambda_tmp1=lambda(count)+(Phi(count)-epsilon*(-1))/(K(count,count));
  lambda_tmp2=lambda(count)+(Phi(count)-epsilon*(1))/(K(count,count));
%         clipping
% xi_tmp=(abs(Phi(count))-epsilon);
% clipper=C*xi_tmp;
%     if(lambda_tmp1<-clipper)
%         lambda_tmp1=-clipper;
%     elseif(lambda_tmp1>clipper)
%         lambda_tmp1=clipper;
%     end
% %        clipping
%     if(lambda_tmp2<-clipper)
%         lambda_tmp2=-clipper;
%     elseif(lambda_tmp2>clipper)
%         lambda_tmp2=clipper;
%     end
    delta1=(lambda_tmp1-lambda(count))*((-1/2)*(lambda_tmp1-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lambda_tmp1)-abs(lambda(count)));
%     -(1/C)*(lambda_tmp1^2-(lambda(count)^2));
    delta2=(lambda_tmp2-lambda(count))*((-1/2)*(lambda_tmp2-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lambda_tmp2)-abs(lambda(count)));
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
K_tmp=K(First,First)*K(Second,Second)-(K(First,Second))^2;
%update lambda(First)
lambda_tmp1=lambda(First)+(Phi(First)*K(Second,Second)-Phi(Second)...
    *K(First,Second)-epsilon*(sgn(First)*K(Second,Second)-sgn(Second)*K(First,Second)))/K_tmp;
% %clipping 
% clipper=C*abs((Phi(First))-epsilon);
%     if(lambda_tmp1<-clipper)
%         lambda_tmp1=-clipper;
%     elseif(lambda_tmp1>clipper)
%         lambda_tmp1=clipper;
%     end
sigma(First)=lambda_tmp1-lambda(First);
sigma1_abs=abs(lambda_tmp1)-abs(lambda(First));
lambda(First)=lambda_tmp1;

%update lambda(Second)
lambda_tmp2=lambda(Second)+(Phi(Second)*K(First,First)-Phi(First)...
    *K(First,Second)-epsilon*(sgn(Second)*K(First,First)...
    -sgn(First)*K(First,Second)))/K_tmp;
% %clipping 
% clipper=C*abs((Phi(Second))-epsilon);
%     if(lambda_tmp2<-clipper)
%         lambda_tmp2=-clipper;
%     elseif(lambda_tmp2>clipper)
%         lambda_tmp2=clipper;
%     end
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
Chi_tmp=abs(Phi_new)-epsilon;
Chi_tmp=Chi_tmp(find(Chi_tmp>0));
Chi(iteration)=norm(Chi_tmp)^(2);
Delta_Chi=Chi(iteration)-Chi(iteration-1);
% NoiseTerm=1/C*NoiseTerm;
%calculate objective function value
% Theta_tmp1=-(1/2)*lambda'*K*lambda+y'*lambda-epsilon*norm(lambda, 1)-0.5*C*Chi(iteration);
Delta_Theta=(-1/2)*(sigma(First)^(2)*K(First, First)+sigma(Second)^(2)*K(Second,Second)+2*sigma(First)*sigma(Second)*K(First,Second))...
    +Phi(First)*sigma(First)+sigma(Second)*Phi(Second)-epsilon*(sigma1_abs+sigma2_abs);
Delta_Theta=Delta_Theta-0.5*C*Delta_Chi;
Theta(iteration)=Theta(iteration-1)+Delta_Theta;
%calculate duality gap
%  G_tmp=y'*lambda-epsilon*(norm(lambda,1))-2*Theta(iteration);
Delta_G=y(First)*sigma(First)+y(Second)*sigma(Second)-epsilon*(sigma1_abs+sigma2_abs)-2*Delta_Theta;
G(iteration)=G(iteration-1)+Delta_G;
if(ratio(iteration-1)<1e-2)
controller=1e-2;
else
    controller=1;
end
ratio(iteration)=controller*abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));
Phi=Phi_new;


end
W_e_real=x'*lambda;
W_e=complex(W_e_real(1:N/2), W_e_real(N/2+1:N));





end

