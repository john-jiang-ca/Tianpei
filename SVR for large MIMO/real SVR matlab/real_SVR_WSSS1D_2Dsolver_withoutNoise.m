function [ symOut, lamida, Theta,G,  iteration, MSE ] = real_SVR_WSSS1D_2Dsolver_withoutNoise( pH,  y, SNRd, symbolContellation, Nr, Nt, M)
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
epsilon=1e-7;   %epsilon
tol=1e-3;  %tolerance of duality gap
K=pH*pH';
C=1;
Phi=zeros(Nr,1);
sigma=zeros(Nr,1);
lamida=zeros(Nr,1);
%% Initialization
Phi=y; %intermediate variable
G=[]; %duality gap
G(1)=0;
for count=1:Nr

    if(abs(Phi(count))>epsilon)
    G(1)=G(1)+(abs(Phi(count))-epsilon);
    end
    
end
G(1)=G(1)*C;
G(1)=1;
Theta(1)=0;
%% update lamida

iteration=1;
while(G(iteration)>tol)

delta=zeros(Nr,1);
sigma_tmp=zeros(Nr, 1);
 %% approximate 2-D work set selection (without damping-bubble sorting)
for count=1:Nr
    lamida_tmp1=lamida(count)+(Phi(count)-epsilon*(-1))/K(count,count);
    lamida_tmp2=lamida(count)+(Phi(count)-epsilon*(1))/K(count,count);
    delta1=(lamida_tmp1-lamida(count))*((-1/2)*(lamida_tmp1-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp1)-abs(lamida(count)));
    delta2=(lamida_tmp2-lamida(count))*((-1/2)*(lamida_tmp2-lamida(count))*K(count,count)+Phi(count))-epsilon*(abs(lamida_tmp2)-abs(lamida(count)));
    if(delta1>delta2)
        delta(count)=delta1;
        lamida_tmp=lamida_tmp1;
        
    else
        delta(count)=delta2;
        lamida_tmp=lamida_tmp2;
    end
    sigma_tmp(count)=lamida_tmp-lamida(count);
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


    


%% update lamida, Phi, Theta and duality gap
sigma(First)=sigma_tmp(First);
sigma(Second)=sigma_tmp(Second);
lamida(First)=lamida(First)+sigma(First);
lamida(Second)=lamida(Second)+sigma(Second);
for count=1:Nr
    Phi(count)=Phi(count)-sigma(First)*K(count,First)-sigma(Second)*K(count,Second);
end

iteration=iteration+1;
%calculate objective function value
Theta(iteration)=-(1/2)*lamida'*K*lamida+y'*lamida-epsilon*norm(lamida, 1);
%calculate duality gap
G(iteration)=y'*lamida-epsilon*(norm(lamida,1))-2*Theta(iteration);
for count=1:Nr
    if(abs(Phi(count))>epsilon)
        G(iteration)=G(iteration)+C*(abs(Phi(count))-epsilon);
    end
end
G(iteration)=G(iteration)/(G(iteration)+Theta(iteration));

end
symOut=(lamida'*pH)';


%% rounding 
d=sqrt(3/(Nt*(M^2-1))); %the distance of constellation
for count=1:Nt
if(M==4)
    if(symOut(count)<=-2*d)
        symOut(count)=-3*d;
    elseif(symOut(count)>-2*d&&symOut(count)<=0)
        symOut(count)=-1*d;
    elseif(symOut(count)>0&&symOut(count)<=2*d)
        symOut(count)=d;
    elseif(symOut(count)>=2*d)
        symOut(count)=3*d;
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
MSE=norm(y-pH*symOut);
end

