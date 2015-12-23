function [ bitOut ] = grayDecoder( symOut, errorV, graycode, symConstell )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
l=length(errorV); % the number of symbol errors
bitOut=zeros(l,1);
for count=1:l
%     for count1=1:length(symConstell)
%         if(abs(symConstell(count1)-symOut(count))<1e-4)
%             index=count1;
%             break
%         end
%     end
    symError=abs(symConstell-symOut(errorV(count)));
    index=find(symError<1e-5);
    bitOut(count)=graycode(index);
end

end