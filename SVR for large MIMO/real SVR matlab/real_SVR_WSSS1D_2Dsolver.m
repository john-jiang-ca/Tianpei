function [ symOut, lamida, Theta,G,  iteration, MSE , NumReliable] = real_SVR_WSSS1D_2Dsolver( pH,  y, SNRd, symbolContellation, Nr, Nt, M,d, epsilon ,C , tol,tau)
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
K=pH*pH';
% C=5;
Phi=zeros(Nr,1);
sigma=zeros(Nr,1);
lamida=zeros(Nr,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
Theta(1)=0;
 for count=1:Nr 
     if(abs(Phi(count))>epsilon)
     Theta(1)=Theta(1)+(abs(Phi(count))-epsilon)^2;
     end
 end
Theta(1)=-0.5*Theta(1)*C;
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
  lamida_tmp1=lamida(count)+(Phi(count)-epsilon*(-1))/(K(count,count));
  lamida_tmp2=lamida(count)+(Phi(count)-epsilon*(1))/(K(count,count));
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
    delta1=(lamida_tmp1-lamida(count))*((-1/2)*(lamida_tmp1-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp1)-abs(lamida(count)));
%     -(1/C)*(lamida_tmp1^2-(lamida(count)^2));
    delta2=(lamida_tmp2-lamida(count))*((-1/2)*(lamida_tmp2-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp2)-abs(lamida(count)));
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
    best_gain=delta(1);
    First=1;
for  count1=2:Nr
    
        if (delta(count1)>best_gain)
            best_gain=delta(count1);
            First=count1;
           
            
        end

end

best_gain=0;
Second=0;
for count1=1:Nr
    if(count1==First)
        continue;
    end
    if(delta(count1)>best_gain)
        best_gain=delta(count1);
        Second=count1;
       
    end
end


    


%% update lamida, sigma, Phi, Theta and duality gap
K_tmp=K(First,First)*K(Second,Second)-(K(First,Second))^2;
%update lamida(First)
lamida_tmp1=lamida(First)+(Phi(First)*K(Second,Second)-Phi(Second)...
    *K(First,Second)-epsilon*(sgn(First)*K(Second,Second)...
    -sgn(Second)*K(First,Second)))/K_tmp;
sigma(First)=lamida_tmp1-lamida(First);
lamida(First)=lamida_tmp1;

%update lamida(Second)
lamida_tmp2=lamida(Second)+(Phi(Second)*K(First,First)-Phi(First)...
    *K(First,Second)-epsilon*(sgn(Second)*K(First,First)...
    -sgn(First)*K(First,Second)))/K_tmp;
sigma(Second)=lamida_tmp2-lamida(Second);
lamida(Second)=lamida_tmp2;
for count=1:Nr
    Phi(count)=Phi(count)-sigma(First)*K(count,First)-sigma(Second)*K(count,Second);
end

iteration=iteration+1;
for count1=1:Nr
     if(abs(Phi(count1))-epsilon>0)
     NoiseTerm=NoiseTerm+(abs(Phi(count1))-epsilon)^2;
     end
% NoiseTerm=NoiseTerm+(lamida(count))^2;
end
% NoiseTerm=1/C*NoiseTerm;
%calculate objective function value
Theta_tmp=-(1/2)*lamida'*K*lamida+y'*lamida-epsilon*norm(lamida, 1);
Theta(iteration)=Theta_tmp-0.5*C*NoiseTerm;
%calculate duality gap
G(iteration)=y'*lamida-epsilon*(norm(lamida,1))-2*Theta(iteration);
G(iteration)=G(iteration)/abs(G(iteration)+Theta(iteration));



end
symOut=(lamida'*pH)';

MSE=norm(y-pH*symOut);
%   symOut_MMSE=zeros(Nt,1);
%   I=ones(Nt,Nt);
%   symOut_MMSE=(inv(pH'*pH)+SNRd^(-1)*I)*pH'*y;
%% rounding 
NumReliable=zeros(Nt,1);
for count=1:Nt
if(M==4)
    if(symOut(count)<=-2*d)
        if(abs(symOut(count)+3*d)<tau*d)
            NumReliable(count)=1;
        end
         symOut(count)=-3*d;
%           symOut_MMSE(count)=-3*d;
    elseif(symOut(count)>-2*d&&symOut(count)<=0)
        if(abs(symOut(count)+d)<tau*d)
            NumReliable(count)=1;
        end
         symOut(count)=-1*d;
%          symOut_MMSE(count)=-1*d;
    elseif(symOut(count)>0&&symOut(count)<=2*d)
                if(abs(symOut(count)-d)>tau*d)
            NumReliable(count)=1;
        end
         symOut(count)=d;
%          symOut_MMSE(count)=d;
    elseif(symOut(count)>=2*d)
          if(abs(symOut(count)-3*d)>tau*d)
            NumReliable(count)=1;
        end
         symOut(count)=3*d;
%           symOut_MMSE(count)=3*d;
    end
elseif (M==8)
        if(symOut(count)<=-6*d)
        symOut(count)=-7*d;
    elseif(symOut(count)>-6*d&&symOut(count)<=-4*d)
        symOut(count)=-5*d;
    elseif(symOut(count)>-4*d&&symOut(count)<=-2*d)
        symOut(count)=-3*d;
    elseif(symOut(count)>-2*d&&symOut(count)<=0)
        symOut(count)=-1*d;
    elseif(symOut(count)>0&&symOut(count)<=2*d)
        symOut(count)=1*d;
    elseif(symOut(count)>2*d&&symOut(count)<=4*d)
        symOut(count)=3*d;
    elseif(symOut(count)>4*d&&symOut(count)<=6*d)
        symOut(count)=5*d;
     elseif(symOut(count)>=6*d)
         symOut(count)=7*d;
        end
    
end
end

%% mean square error
% MSE=norm(y-pH*symOut);
end

