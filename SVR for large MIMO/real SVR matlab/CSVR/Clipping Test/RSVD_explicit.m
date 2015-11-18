function [ symOut] = RSVD_explicit( H,  y, SNRd,  M, pav, C, tol , epsilon)
%real SVR detecor  
%OUTPUT
% symOut: the estimated regression coefficience vector (sliced by symbol constellation alphabet) 


%Input
% H: data sample matrix (input) 
% y: data sample vector (output)
% M: M-QAM 
% pav: average power of the transmit symbols
% SNRd: receive signal to noise ratio 
% C: constant to control the tradeoff between regularization term and
% outlier penaty term
% epsilon: epsilon for precision control 
% tol: tolerance of the duality gap that determine when the algorithm stops  




N=length(H(1,:));  %the dimentsion of input 
L=length(H(:,1));  %the number of input
lambda=zeros(L,1);  %dual optimization variable vector
sigma=zeros(L,1); %sigma=lambda_new-lambda
Phi=zeros(L,1);   %intermediate variable Phi
G=[]; %duality gap
Theta=[]; %dual objective function
Chi=[];  %noise Term
ratio=[];  %ratio=abs(G)/(abs(G)+abs(Theta));
iteration=1; %iteration time
%% Initialization
K=H*H';   %calculate kernle matrix
Phi=y; %intermediate variable
G(1)=0;
Theta(1)=0;
Phi_tmp=abs(Phi)-epsilon;
Phi_tmp=Phi_tmp(find(Phi_tmp>0));
Chi(1)=norm(Phi_tmp)^(2);
Theta(1)=-0.5*C*Chi(1);
G(1)=-2*Theta(1);
ratio(1)=G(1)/(G(1)+Theta(1));
%% Quadratic Programming of the dual objective function 
while(ratio(iteration)>tol)   %performs update to dual optimization variables lambda until the duality gap satisfies the tolerance

delta=zeros(L,1);  %gain of the sub dual objective function 
sgn=zeros(L,1);  %record the sign of update value of dual optimizaition variables
%% single direction searching 
for count=1:L
  lambda_tmp1=lambda(count)+(Phi(count)-epsilon*(-1))/(K(count,count));
  lambda_tmp2=lambda(count)+(Phi(count)-epsilon*(1))/(K(count,count));
    %% explicit clipping
    sigma_tmp1=lambda_tmp1-lambda(count);
    sigma_tmp2=lambda_tmp2-lambda(count);
    Phi_tmp1=Phi(count)-K(count,count)*sigma_tmp1;
% Phi_tmp1=Phi(count);
% Phi_tmp2=Phi(count);
    Phi_tmp2=Phi(count)-K(count,count)*sigma_tmp2;
%   bottom_tmp1=-C*max(0, -Phi_tmp1);
%   upper_tmp1=C*max(0, Phi_tmp1);
%   bottom_tmp2=-C*max(0, -Phi_tmp2);
%   upper_tmp2=C*max(0, Phi_tmp2);
%   if(lambda_tmp1<bottom_tmp1)
%       lambda_tmp1=bottom_tmp1;
%   elseif (lambda_tmp1>upper_tmp1)
%       lambda_tmp1=upper_tmp1;
%   end
%   
%   if(lambda_tmp2<bottom_tmp2)
%       lambda_tmp2=bottom_tmp2;
%   elseif (lambda_tmp2>upper_tmp2)
%       lambda_tmp2=upper_tmp2;
%   end
  
  delta1=(lambda_tmp1-lambda(count))*((-1/2)*(lambda_tmp1-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lambda_tmp1)-abs(lambda(count)));
  delta2=(lambda_tmp2-lambda(count))*((-1/2)*(lambda_tmp2-lambda(count))*K(count,count)+Phi(count))-epsilon*(abs(lambda_tmp2)-abs(lambda(count)));

  
  
    if(delta1>delta2)
        delta(count)=delta1;
        sgn(count)=-1;
        
    else
        delta(count)=delta2;
        sgn(count)=1;
    end

end
%find the coordinates whose updates can generate first and second largest
%gain to sub dual objective function (indexed by 'First' and 'Second')
[value, First]=max(delta);
delta(First)=min(delta)-1;
[value, Second]=max(delta);


    



%% double Direction update
%Performs double Direction update of the 'First' and 'Second' dual
%optimization variables that determined in Single direction searching
%process
K_tmp=K(First,First)*K(Second,Second)-(K(First,Second))^2;

lambda_tmp1=lambda(First)+(Phi(First)*K(Second,Second)-Phi(Second)...
    *K(First,Second)-epsilon*(sgn(First)*K(Second,Second)-sgn(Second)*K(First,Second)))/K_tmp;

lambda_tmp2=lambda(Second)+(Phi(Second)*K(First,First)-Phi(First)...
    *K(First,Second)-epsilon*(sgn(Second)*K(First,First)...
    -sgn(First)*K(First,Second)))/K_tmp;

%% clipping
    sigma_tmp1=lambda_tmp1-lambda(First);
    sigma_tmp2=lambda_tmp2-lambda(Second);
    Phi_tmp1=Phi(First)-K(First,First)*sigma_tmp1-K(First, Second)*sigma_tmp2;
    Phi_tmp2=Phi(Second)-K(Second, First)*sigma_tmp1-K(Second, Second)*sigma_tmp2;
% Phi_tmp1=Phi(First);
% Phi_tmp2=Phi(Second);
%   bottom_tmp1=-C*max(0, -Phi_tmp1);
%   upper_tmp1=C*max(0, Phi_tmp1);
%   bottom_tmp2=-C*max(0, -Phi_tmp2);
%   upper_tmp2=C*max(0, Phi_tmp2);
%   if(lambda_tmp1<bottom_tmp1)
%       lambda_tmp1=bottom_tmp1;
%   elseif (lambda_tmp1>upper_tmp1)
%       lambda_tmp1=upper_tmp1;
%   end
%   
%   if(lambda_tmp2<bottom_tmp2)
%       lambda_tmp2=bottom_tmp2;
%   elseif (lambda_tmp2>upper_tmp2)
%       lambda_tmp2=upper_tmp2;
%   end


%update lambda(First)
sigma(First)=lambda_tmp1-lambda(First);
sigma1_abs=abs(lambda_tmp1)-abs(lambda(First));
lambda(First)=lambda_tmp1;
%update lambda(Second)
sigma(Second)=lambda_tmp2-lambda(Second);
sigma2_abs=abs(lambda_tmp2)-abs(lambda(Second));
lambda(Second)=lambda_tmp2;

iteration=iteration+1;
%% update Chi, Phi, Theta, G and ratio
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
ratio(iteration)=controller*abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));  %update ratio
Phi=Phi_new;  %update Phi


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