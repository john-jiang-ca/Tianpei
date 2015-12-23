function [ H1, H2, Index1, Index2 ] = Antenna_Sel(H, N, SNR )
%this function implement channel selection based on diversity maximization
%antenna selection scheme.
%INPUT
%H: input channel matrix
%N: number of channel selected for the first stage
%SNR: symbol signal to noise ratio
%OUTPUT
%H1: n columns detected using brute-force searching
%H2: Nt-n column detected using MMSE/MMSE-OSIC
%Index1: antennas indexes chosen in the first stage
%Index2: antennas indexes chosen in the second stage
%Tianpei Chen
%Dec 02 2015
Nr=length(H(:,1));   % the number of receive antennas
Nt=length(H(1,:));   %the number of transmit antennas
candi=combntns([1:Nt],N);      %choose N antennas outof Nt
size=length(candi(:,1));      %the size of the candidate list
I=eye(Nt-N);    %the identity matrix
maxD=zeros(size,1);  %used to store all the max diagonal value of the inverse matrix
for count=1:size
    H2_tmp=H;
    H2_tmp(:, candi(count,:))=[];
    G=I/(H2_tmp'*H2_tmp+(SNR^(-1))*I);
    [maxD(count), nonU]=max(diag(G));
    
end
[value, index]=min(maxD);    %find the index of the candidate that can be chosen
Index1=candi(index,:);     % the index vector of the channel selected in the  first stage
Index2=1:Nt;
Index2(:,Index1)=[];   %the index vector of the channel selected in the second stage
H1=H(:, Index1);   %construct H1
H2=H(:, Index2);   %construct H2

    


end