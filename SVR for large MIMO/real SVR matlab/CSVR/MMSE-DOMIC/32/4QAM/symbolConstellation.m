function [ symConstell ] = symbolConstellation( M, pav )
%Generage symbol constellation
%   Detailed explanation goes here
Q=sqrt(M);
symConstell=zeros(M,1);
symConstellReal=zeros(Q,1);
symReal=zeros(Q,1);
symImag=zeros(Q,1);
if(M==2)
  d=sqrt(pav);  
else
    d=sqrt(3*pav/(2*(M-1)));
end

for count=1:Q
    symConstellReal(count)=-(Q-1)*d+(count-1)*2*d;
end
for count=1:M
    symReal(count)=symConstellReal(fix((count-1)/Q)+1);
    symImag(count)=symConstellReal(Q-mod(count-1,Q));
    symConstell(count)=complex(symReal(count),symImag(count));
end

