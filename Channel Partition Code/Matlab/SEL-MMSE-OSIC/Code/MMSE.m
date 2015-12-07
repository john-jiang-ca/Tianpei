function [ symOut, outlier ] = MMSE( y, H, SNRd, M, pav )
%Minimum Mean Square Error Detector
% Rectangular_QAM_slicer performs slicing based on the symbol constellation
% alphabet
Nt=length(H(1,:));
I=eye(Nt);
G=(I/(H'*H+(SNRd/Nt)^(-1)*I))*H';
symOut_tmp=G*y;
symOut=Rectangular_QAM_slicer(symOut_tmp, M, pav);
end

