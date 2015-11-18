function [ symOut, outlier ] = MMSE( y, H, SNRd, M, pav )
%Minimum Mean Square Error Detector
% Rectangular_QAM_slicer performs slicing based on the symbol constellation
% alphabet
Nt=length(H(1,:));
I=eye(Nt);
G=(I/(H'*H+SNRd^(-1)*I))*H';
symOut_tmp=G*y;
if(M==2)
    d=sqrt(pav);
else
d=sqrt(3*pav/(2*(M-1)));
end
Q=ceil(sqrt(M));
symOut_tmp_real=[real(symOut_tmp); imag(symOut_tmp)];
outlier=length(find(symOut_tmp_real>(Q-1)*d))+length(find(symOut_tmp_real<-(Q-1)*d));
symOut=Rectangular_QAM_slicer(symOut_tmp, M, pav);
end

