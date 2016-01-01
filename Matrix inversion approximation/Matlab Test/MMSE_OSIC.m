function [ symOut ] = MMSE_OSIC(y, H, SNR, M, pav )
% This subroutine performs MMSE_OSIC algorithm
%INPUT:
%y: received signal vector
%H: channel propagation matrix
%SNR: symbol signal to noise ratio
%OUTPUT:
%symOut: estimation of the transmit symbol vector

Nr=length(H(:,1));
Nt=length(H(1,:));
sorting=1:Nt; %location list that store the original location of the detected symbol
H_tmp=H; %temporary channel matrix
y_tmp=y; %residual observation of SIC
symOut=zeros(Nt,1);  %the output estimation symbol vector
for count=1:Nt
    I=eye(length(H_tmp(1,:)));   %identity matrix 
    G_pre=I/(H_tmp'*H_tmp+(SNR^(-1))*I);  
    [value, index]=min(diag(G_pre));  %choose the symbol with the strongest postprocessing SNR 
    s_tmp=G_pre(index, :)*H_tmp'*y_tmp;   %get the estimation of the symbol chosen
    s_tmp=Rectangular_QAM_slicer(s_tmp, M, pav); %slicing based on the signla constellation
    symOut(sorting(index))=s_tmp;  %put the detected symbol to the original location of the transmit symbol vector
    y_tmp=y_tmp-H_tmp(:,index)*s_tmp; %perform interference cancellation
    H_tmp(:,index)=[]; %update channel matrix
    sorting(index)=[]; %update location list
end




end