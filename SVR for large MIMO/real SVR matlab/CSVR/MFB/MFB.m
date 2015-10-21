function [ symOut ] = MFB( y, H, symConstell )
%match filter bound
%   Detailed explanation goes here
M=length(symConstell);
Euclidean=zeros(M,1);
for count=1:M
Euclidean(count)=norm(y-H*symConstell(count));
end
[value, index]=min(Euclidean);
symOut=symConstell(index);

end

