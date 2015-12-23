function [ value ] = pdfSub( Nr, Nt, rho, candidate )
% This subroutine calculate the sub function of pdf of orthogonality
% measure
% 
Cvalue=1;
for i=2:Nt
    Cvalue=Cvalue*((-1)^(candidate(i-1))*nchoosek((i-2), candidate(i-1))...
        /((i-1+candidate(i-1))*beta((Nr-i+1), (i-1))));
end

end

