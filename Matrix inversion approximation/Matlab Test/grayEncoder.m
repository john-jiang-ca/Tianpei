function [ graycode ] = grayEncoder( M )
%gray encoder 
%M is the constellation size
if(M==4)
    graycode=[0,1,2,3];
elseif(M==16)
    graycode=[0,1,3,2,4,5,7,6,12,13,15,14,8,9,11,10];
end


end