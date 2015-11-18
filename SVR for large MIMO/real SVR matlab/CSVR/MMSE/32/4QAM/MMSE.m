function [ symOut, MSE ] = MMSE( y, H, SNRd, M, pav, dataMod )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
Nr=length(H(:,1));
Nt=length(H(1,:));
I=eye(Nt);
G=(I/(H'*H+SNRd^(-1)*I))*H';
symOut_tmp=G*y;
MSE=norm(dataMod-symOut_tmp)^(2)/Nr;
symOut=Rectangular_QAM_slicer(symOut_tmp, M, pav);

end

