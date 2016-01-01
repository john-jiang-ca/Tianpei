function [ symOut ] = MMSE_A( y, H, SNRd, M, pav )
%Minimum Mean Square Error Detector
% Rectangular_QAM_slicer performs slicing based on the symbol constellation
% alphabet
Nt=length(H(1,:));
I=eye(Nt);
% G=(I/(H'*H+(SNRd)^(-1)*I))*H';
G_inter=H'*H+SNRd^(-1)*I;
[G_interI, nonOfuse]=SENIA(G_inter);
G=G_interI*H';
symOut_tmp=G*y;
symOut=Rectangular_QAM_slicer(symOut_tmp, M, pav);
end