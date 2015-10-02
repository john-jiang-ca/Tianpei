function [  symOut,   MSE] = MMSE(pH,  y, SNRd, symbolContellation, Nr, Nt, M )
%MMSE algorithm
%Output: symOut, and mean square error
G=zeros(Nt,Nr);
symOut=zeros(Nt,1);
I=ones(Nt,Nt);
G=(inv(pH'*pH+SNRd^(-1)*I))*pH';
symOut=G*y;
MSE=norm(y-pH*symOut);
%% rounding 
d=sqrt(3/(Nt*(M^2-1))); %the distance of constellation
for count=1:Nt
if(M==4)
    if(symOut(count)<=-2*d)
         symOut(count)=-3*d;
%           symOut_MMSE(count)=-3*d;
    elseif(symOut(count)>-2*d&&symOut(count)<=0)
         symOut(count)=-1*d;
%          symOut_MMSE(count)=-1*d;
    elseif(symOut(count)>0&&symOut(count)<=2*d)
         symOut(count)=d;
%          symOut_MMSE(count)=d;
    elseif(symOut(count)>=2*d)
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
%% Mean Square Error
% MSE=norm(y-pH*symOut);

end

