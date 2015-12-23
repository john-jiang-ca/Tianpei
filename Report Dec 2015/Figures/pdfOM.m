function [ value ] = pdfOM( Nr ,Nt )
% This routine calculate the pdf of orthogonality measure
% INPUT: 
%Nr: number of receive antennas
%Nt: number of transmit antennas
%OUTPUT
%value: the value of the pdf
%% parameters settings
step=0.01;
numberT=round(1/step)+1;   %the total number of points considered
value=zeros(numberT);    %pdf value vector
candidateV=1:1:(Nt-1);   
candidate=fullfact(candidateV); %% generate all the possible ji combinations
%% sub function





end

