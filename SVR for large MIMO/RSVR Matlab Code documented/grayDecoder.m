function [ bitOut ] = grayDecoder( symOut, graycode, symConstell )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
Nt=length(symOut);
bitOut=zeros(Nt,1);
for count=1:Nt
    for count1=1:length(symConstell)
        if(abs(symConstell(count1)-symOut(count))<1e-4)
            index=count1;
            break
        end
    end
    bitOut(count)=graycode(index);
end

end

