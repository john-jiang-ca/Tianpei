function [index1,index2, sigma1, sigma2, sigma1_a, sigma2_a, Lambda_new, Phi_new, deltaTheta]=a2DSolver_modified(Kr, Lambda, Phi, epsilon)
%this subroutine perform sequential  or on shot update of the optimization
%variables
Nr=length(Kr(1,:));
%% perform one shot searching
% sigma_single_v=zeros(Nr,1);
% sigma_a_single_v=zeros(Nr,1);
deltaTheta_v=zeros(Nr,1);
% Lambda_new_v=zeros(Nr,1);
sgn=zeros(Nr,1);
for count=1:Nr
    lambda_tmp1=Lambda(count)+(Phi(count)-epsilon*1)/2*Kr(count,count);   %different from the original
    lambda_tmp2=Lambda(count)+(Phi(count)-epsilon*(-1))/2*Kr(count,count);
    sigma_tmp1=lambda_tmp1-Lambda(count);
    sigma_tmp1_a=abs(lambda_tmp1)-abs(Lambda(count));
    sigma_tmp2=lambda_tmp2-Lambda(count);
    sigma_tmp2_a=abs(lambda_tmp2)-Lambda(count);
    deltaThetaSingle_tmp1=-((sigma_tmp1)^(2))*Kr(count,count)+Phi(count)*(sigma_tmp1)-epsilon*(sigma_tmp1_a); %different from the original
    deltaThetaSingle_tmp2=-((sigma_tmp2)^(2))*Kr(count,count)+Phi(count)*(sigma_tmp2)-epsilon*(sigma_tmp2_a);
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
Lambda_new=Lambda;

%% perform sequential searching
%% perform 2D update
danomin_tmp=2*(Kr(index1,index1)*Kr(index2,index2)-(Kr(index1,index2)^(2)));  %different from the original
Lambda_new1=Lambda(index1)+(Phi(index1)*Kr(index2,index2)-Phi(index2)*Kr(index1,index2)...
    -epsilon*(sgn1*Kr(index2,index2)-sgn2*Kr(index1,index2)))/danomin_tmp;
Lambda_new2=Lambda(index2)+(Phi(index2)*Kr(index1,index1)-Phi(index1)*Kr(index1,index2)...
    -epsilon*(sgn2*Kr(index1,index1)-sgn1*Kr(index1,index2)))/danomin_tmp;
sigma1=Lambda_new1-Lambda(index1);
sigma2=Lambda_new2-Lambda(index2);
sigma1_a=abs(Lambda_new1)-abs(Lambda(index1));
sigma2_a=abs(Lambda_new2)-abs(Lambda(index2));
Lambda_new(index1)=Lambda_new1;
Lambda_new(index2)=Lambda_new2;
%% update Lambda, Phi, Psi, and deltaTheta

% sigma1=sigma_single_v(index1);
% sigma2=sigma_single_v(index2);
% sigma1_a=sigma_a_single_v(index1);
% sigma2_a=sigma_a_single_v(index2);
% Lambda_new(index1)=Lambda(index1)+sigma1;
% Lambda_new(index2)=Lambda(index2)+sigma2;
% deltaTheta=deltaTheta1+deltaTheta2-sigma1*sigma2*Kr(index1,index2);
deltaTheta=-(sigma1^(2)*Kr(index1,index1)+sigma2^(2)*Kr(index2,index2)+2*sigma1*sigma2*Kr(index1,index2))+Phi(index1)*sigma1+Phi(index2)*sigma2...
    -epsilon*(sigma1_a+sigma2_a);
Phi_new=Phi-2*sigma1*Kr(:,index1)-2*sigma2*Kr(:,index2);
% if(label==1)
%     Psi_new=Psi-sigma1*Ki(:,index1)-sigma2*Ki(:,index2);
% else
%     Psi_new=Psi+sigma1*Ki(:,index1)+sigma2*Ki(:,index2);
% end

end