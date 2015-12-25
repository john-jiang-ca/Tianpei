function [ value ] = pdfSub( k1, k2, candidate, rho )
% This subroutine calculate the sub function of pdf of orthogonality
% measure
% 
Cvalue=1;
L=length(k1);  %Nt-1
for count=1:L
    Cvalue=Cvalue*((-1)^(candidate(count))*nchoosek(k2(count)-1, candidate(count))...
        /(beta(k1(count), k2(count))));
end
value=0;
for count1=1:L
    prodTmp=1;
    for count2=1:L
        if(count2==count1)
            continue;
        end
        prodTmp=prodTmp*(k1(count2)+candidate(count2)-k1(count1)-candidate(count1));
    end
    value=value+rho^(k1(count1)+candidate(count1)-1)/prodTmp;
end

end

