function [ value ] = pdfOM( Nr ,Nt )
% This routine calculate the pdf of orthogonality measure
% INPUT: 
%Nr: number of receive antennas
%Nt: number of transmit antennas
%OUTPUT
%value: the value of the pdf
%% parameters settings
step=0.01;
numberT=round(1/step);   %the total number of points considered
value=zeros(numberT,1);    %pdf value vector
candidateV=1:1:(Nt-1);   
candidate=fullfact(candidateV); %% generate all the possible ji combinations
candidate=candidate-ones(size(candidate));
k1=zeros(Nt-1,1);   %shape parameter k1 
k2=zeros(Nt-1,1);   %shape parameter k2
for count1=2:Nt
    k1(count1-1)=Nr-count1+1;
    k2(count1-1)=count1-1;
end
%% sub function
numberCandidate=length(candidate(:,1));   %the number of summation terms
for count1=1:numberT    % the number of sample point of pdf
    for count2=1:numberCandidate   % the number of summation terms
        valueTmp=pdfSub(k1,k2, candidate(count2,:), count1*step);
        value(count1)=value(count1)+valueTmp;
    end
end



end

