function [ bitError ] = checkBitError( bitOut, errorV, grayData,M )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
bitError=0;
l=length(errorV);
for count=1:l
    bit=bitOut(count);
    data=grayData(errorV(count));
    for count1=1:ceil(log2(M))
        bitCurrent=mod(bit, 2);
        dataCurrent=mod(data,2);
        bit=floor(bit/2);
        data=floor(data/2);
       if (bitCurrent~=dataCurrent)
           bitError=bitError+1;
       end
    end
end

end

