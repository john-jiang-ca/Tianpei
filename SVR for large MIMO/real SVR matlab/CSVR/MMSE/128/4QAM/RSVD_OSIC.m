function [ symOut ] = RSVD_OSIC(y, H, SNRd, M, pav)
%this routine implement MMSE-OSIC algorithm
%   Detailed explanation goes here

Nr=length(H);
Nt=length(H(:,1));
symOut=zeros(Nt,1);  %output symbol vector
% symOut_tmp_v=zeros(Nt,1);  %
symOut_list=1:Nt; %record the current index of the updated subvector
% symOut_Recons=zeros(Nt,1);

H_tmp=H;  %updated channel matrix
y_tmp=y; %updated received symbol vector 
for count=1:Nt
    I=eye(Nt-(count-1));
    G=(I/(H_tmp'*H_tmp+SNRd^(-1)*I));   %inverse of HHH matrix with SNR
    [value, index]=min(diag(G));   %choice of the channel with the minimum post processing SNR
%     G=G*H_tmp';   %equalization matrix
%     symOut_tmp=G(index,:)*y_tmp;   
%     symOut_tmp=Rectangular_QAM_slicer(symOut_tmp,M,pav); %detected symbol 
y_tmp_r=[real(y_tmp); imag(y_tmp)];
H_tmp_r=[real(H_tmp), -imag(H_tmp); imag(H_tmp), real(H_tmp)];
[symOut_tmp]=RSVR(H_tmp_r, y_tmp_r, SNRd, M, pav);
    symOut(symOut_list(index))=symOut_tmp(index); %upate output symbol vector
    y_tmp=y_tmp-H_tmp(:,index)*symOut_tmp(index); %update received symbol vector
    H_tmp(:,index)=[];  %update channel matrix
    symOut_list(index)=[];   %update symbol index list

end



end

