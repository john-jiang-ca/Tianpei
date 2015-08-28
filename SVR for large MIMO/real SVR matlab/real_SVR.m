function [ symOut, lamida, Theta,G,  iteration, Error ] = real_SVR( pH,  y, SNRd, symbolContellation, Nr, Nt, M)
%real SVR detection 
%   Output
% symOut: regression coefficiences estimation
% Theta: the vector that store all the objective function value
%Input
% pH: data sample
% y: data  sample output
% SNRd: signal to noise ratio at output 
% symbolConstellation: the symbol constellation
% iteration: iteration time
epsilon=1e-2;   %epsilon
tol=1e-2;  %tolerance of duality gap
K=pH*pH';
Phi=zeros(Nr,1);
lamida=zeros(Nr,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
for count=1:Nr

    if(abs(Phi(count))>epsilon)
    G(1)=G(1)+(abs(Phi(count))-epsilon)^2;
    end
end
Theta(1)=0;
%% update lamida
First=-1;
Second=First;
bestGain=0;
delta=0;
lamida1=0;
lamida2=0;
sigma1=0;
sigma2=0;
iteration=1;
while(G(iteration)>tol)
    
    
 %% approximate 2-D work set selection (without damping)
for count=1:Nr
    lamida_tmp1=lamida(count)+(Phi(count)-epsilon*(-1))/K(count,count);
    lamida_tmp2=lamida(count)+(Phi(count)-epsilon*(1))/K(count,count);
    delta1=(lamida_tmp1-lamida(count))*((-1/2)*(lamida_tmp1-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp1)-abs(lamida(count)));
    delta2=(lamida_tmp2-lamida(count))*((-1/2)*(lamida_tmp2-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp2)-abs(lamida(count)));
    if(delta1>delta2)
        delta=delta1;
        lamida_tmp=lamida_tmp1;
        
    else
        delta=delta2;
        lamida_tmp=lamida_tmp2;
    end
    sigma_tmp=lamida_tmp-lamida(count);
    if(delta>bestGain)
        bestGain=delta;
        Second=First;
        lamida2=lamida1;
        sigma2=sigma1;
        First=count;
        lamida1=lamida_tmp;
        sigma1=sigma_tmp;
    end
    
end

%% update lamida, Phi, Theta and duality gap
lamida(First)=lamida1;
lamida(Second)=lamida2;
for count=1:Nr
    Phi(count)=Phi(count)-sigma1*K(count,First)-sigma2*K(count,Second);
end

iteration=iteration+1;
%calculate objective function value
Theta(iteration)=-(1/2)*lamida'*K*lamida+y'*lamida-epsilon*norm(lamida, 1);
%calculate duality gap
G(iteration)=y'*lamida-epsilon*(norm(lamida,1))-2*Theta(iteration);
for count=1:Nr
    if(abs(Phi(count))>epsilon)
        G(iteration)=G(iteration)+(abs(Phi(count))-epsilon)^2;
    end
end


end
symOut=lamida'*pH;
%% mean square error
Error=norm(y-pH*symOut);
%% rounding 
d=sqrt(3/(Nt*(M^2-1))); %the distance of constellation
for count=1:Nt
if(M==4)
    if(symOut(count<=-2*d))
        symOut(count)=-3*d;
    elseif(symOut(count)>-2*d&&symOut(count)<=0)
        symOut(count)=-1*d;
    elseif(symOut(count)>0&&symOut(count)<=2*d)
        symOut(count)=d;
    elseif(symOut(count)>=2*d)
        symOut(count)=3*d;
    end
elseif (M==8)
        if(symOut(count<=-6*d))
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

end

