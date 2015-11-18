function [Lambda]=a2DSolver_modified(y, K, C, epsilon, tol)
%this subroutine perform sequential  or on shot update of the optimization
%variables
%% Initialization
Nr=length(K(1,:));
Lambda=zeros(Nr,1);
Phi=y;
Chi_v=abs(Phi)-epsilon;
Chi_v=Chi_v(find(Chi_v>0));
Chi=norm(Chi_v)^(2);
Theta(1)=-(C/2)*Chi;
G(1)=-2*Theta;
ratio(1)=G(1)/(G(1)+Theta(1));
iteration=1;
while(ratio(iteration)>=tol)
%% perform one shot searching
deltaTheta_v=zeros(Nr,1);
sgn=zeros(Nr,1);
index1=0;
index2=0;
for count=1:Nr
    danomin_tmp1=2*K(count,count);
    lambda_tmp1=Lambda(count)+(Phi(count)-epsilon*1)/danomin_tmp1;   %different from the original
    lambda_tmp2=Lambda(count)+(Phi(count)-epsilon*(-1))/danomin_tmp1;
    sigma_tmp1=lambda_tmp1-Lambda(count);
    sigma_tmp1_a=abs(lambda_tmp1)-abs(Lambda(count));
    sigma_tmp2=lambda_tmp2-Lambda(count);
    sigma_tmp2_a=abs(lambda_tmp2)-Lambda(count);
    deltaThetaSingle_tmp1=-((sigma_tmp1)^(2))*K(count,count)+Phi(count)*(sigma_tmp1)-epsilon*(sigma_tmp1_a); %different from the original
    deltaThetaSingle_tmp2=-((sigma_tmp2)^(2))*K(count,count)+Phi(count)*(sigma_tmp2)-epsilon*(sigma_tmp2_a);
    [value, index]=max([deltaThetaSingle_tmp1, deltaThetaSingle_tmp2]);
    if(index==1)
%         Lambda_new_v(count)=lambda_tmp1;
        sgn(count)=1;
%         sigma_single_v(count)=sigma_tmp1;
%         sigma_a_single_v(count)=sigma_tmp1_a;
        deltaTheta_v(count)=value;
    else
%         Lambda_new_v(count)=lambda_tmp2;
       sgn(count)=-1;
%         sigma_single_v(count)=sigma_tmp2;
%         sigma_a_single_v(count)=sigma_tmp2_a;
        deltaTheta_v(count)=value;  
    end
    
end

[deltaTheta1, index1]=max(deltaTheta_v);
deltaTheta_v(index1)=min(deltaTheta_v)-1;
[deltaTheta2, index2]=max(deltaTheta_v);
sgn1=sgn(index1);
sgn2=sgn(index2);
% Lambda_new=Lambda;

%% perform sequential searching
%% perform 2D update
danomin_tmp=2*(K(index1,index1)*K(index2,index2)-(K(index1,index2)^(2)));  %different from the original
Lambda_new1=Lambda(index1)+(Phi(index1)*K(index2,index2)-Phi(index2)*K(index1,index2)...
    -epsilon*(sgn1*K(index2,index2)-sgn2*K(index1,index2)))/danomin_tmp;
Lambda_new2=Lambda(index2)+(Phi(index2)*K(index1,index1)-Phi(index1)*K(index1,index2)...
    -epsilon*(sgn2*K(index1,index1)-sgn1*K(index1,index2)))/danomin_tmp;
sigma1=Lambda_new1-Lambda(index1);
sigma2=Lambda_new2-Lambda(index2);
sigma1_a=abs(Lambda_new1)-abs(Lambda(index1));
sigma2_a=abs(Lambda_new2)-abs(Lambda(index2));
%calculate the change of Theta_S
deltaTheta_S=-(sigma1^(2)*K(index1,index1)+sigma2^(2)*K(index2,index2)+2*sigma1*sigma2*K(index1,index2))+Phi(index1)*sigma1+Phi(index2)*sigma2...
    -epsilon*(sigma1_a+sigma2_a);
%update Lambda
Lambda(index1)=Lambda_new1;
Lambda(index2)=Lambda_new2;
%% update Lambda, Phi, and Theta, G and ratio
Phi=Phi-2*sigma1*K(:,index1)-2*sigma2*K(:,index2);  %update Phi
Chi_new_v=abs(Phi)-epsilon;
Chi_new_v=Chi_new_v(find(Chi_new_v>0));
Chi_new=norm(Chi_new_v)^(2);
deltaNoise=Chi_new-Chi; %record the change of noise term
Chi=Chi_new;   %update  noise term
deltaTheta=deltaTheta_S-(C/2)*deltaNoise;
deltaG=y(index1)*sigma1+y(index2)*sigma2-epsilon*(sigma1_a+sigma2_a)-2*deltaTheta;
iteration=iteration+1;


Theta(iteration)=Theta(iteration-1)+deltaTheta; %update Theta
% Theta_tmp=-Lambda'*K*Lambda+y'*Lambda-epsilon*(norm(Lambda,1))-(C/2)*Chi;
G(iteration)=G(iteration-1)+deltaG;   %update G
% G_tmp=y'*Lambda-epsilon*norm(Lambda,1)-2*Theta_tmp;
if(ratio(iteration-1)<1e-2)
    controller=1e-2;
else
    controller=1;
end
% controller=1;
ratio(iteration)=controller*abs(G(iteration))/(abs(G(iteration))+abs(Theta(iteration)));
% ratio(iteration)=G(iteration)/(G(iteration)+Theta(iteration));
% if(label==1)
%     Psi_new=Psi-sigma1*Ki(:,index1)-sigma2*Ki(:,index2);
% else
%     Psi_new=Psi+sigma1*Ki(:,index1)+sigma2*Ki(:,index2);
% end
end
end