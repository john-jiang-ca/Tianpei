function [ bitError ] = checkBitError( bitOut, grayData,M )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
bitError=0;
Nt=length(bitOut);
for count=1:Nt
    bit=bitOut(count);
    data=grayData(count);
    for count1=1:ceil(log2(M))
       currentBit=mod(bit ,2);
       currentData=mod(data,2);
       bit=floor(bit/2);
       data=floor(data/2);
       if (currentBit~=currentData)
           bitError=bitError+1;
       end
    end
end

end

