
%% Main System File for Bit Error Rate Monte-Carlo Simulation 
%=============
%Author: Tianpei Chen
% Introduction:
% this routine performs the Monte-Carlo Simulation of the bit error rate of
% proposed detectors in Large Scale-MIMO (LS-MIMO) systems. A set of signal
% to noise ratio (SNR) points are considered, (Notice: the SNR in Monte-Carlo simulation is defined
% as the average receive SNR generally, however here we use the difinition in Dejelili's thesis that is, SNR of each symbol). 
% For each SNR point, the
% simulations works in the step of channel realizations, symbol error rate
% and bit error rate of the detectors are calculated when both a certain number
% of symbol errors are accumulated and adequate channel realization times
% is reached.
% this routine output the BER or SER or related simulation results to txt
% files.
%Date: Dec 02 20015
%=============
%% Function Description
% grayEncoder: generate gray code list
% symbolConstellation: generate symbol constellation alphabet 
% sel_MMSE_OSIC: selection MMSE OSIC detector
% grayDecoder: decoding symbol estimations to corresponding gray code
% checkBitError: comparing the output bit sequence and the original
% transmit bit sequence, then calculate bit errors

close all
clear all
tic

%% System Configuration
M =16;         %size of the signal constellation alphabet  (rectangular M-QAM)
Nt=4;         %number of transmit antennas
Nr=4;         %number of receive antennas
SNR=[0:2:16];       %average receive signal to noie ratio (SNR) in (dB) 
SNRd=10.^(SNR.*0.1);   %SNR in dicimal 
pav=1/Nt;  %average power of the transmitted symbols
% noiseV=pav./SNRd;   %noise variance of  Additive White Gaussian Noise (AWGN)
noiseV=1./SNRd;   
SER_sel_MMSE_OSIC=zeros(length(SNR),1);    %symbol error rate of MMSE detector for every SNR (dB) 
BER_sel_MMSE_OSIC=zeros(length(SNR),1);    %bit error rate of MMSE detector for every SNR (dB)
ChannelRealization=5e3;   %the minimum number of channel realization
minSymbolError=200;       %the minimum number of symbol error accumulated 
%% Signal Modulation
graycode=grayEncoder(M); %gray code encoder
[symConstell]=symbolConstellation( M, pav );  %generate symbol constellation

%% generate file to record the simulation results
fid=fopen('F:\GitHub\Tianpei\Channel Partition Code\Matlab\SEL-MMSE-OSIC\data\MMSE_OSIC.txt', 'a');
fprintf(fid, '\n');
fprintf(fid, '=========================================================================================\n');
fprintf(fid ,'this file compare the BER simulation results of MMSE in LS(large scale)-MIMO systems\n');
fprintf(fid, 'SYSTEM CONFIGURATION\n');
fprintf(fid, '****************\n');
fprintf(fid, '%d X %d MIMO system with %d QAM modulation\n', Nr,Nt,M);
fprintf(fid, 'the average power of transmit symbols is: %0.10f\n', 1/Nt);
fprintf(fid,'the average receive SNR (dB) are:\n');
for count=1:length(SNR)
    fprintf(fid, '%d ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'The minimum channel realization is %e\n', ChannelRealization);
fprintf(fid, 'The minimum symbol error accumulatted is %e\n', minSymbolError);
fprintf(fid, '****************\n');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'SIMULATION OUTPUT\n');
fprintf(fid, '****************\n');
fprintf(fid, 'the first column shows the current SNR (dB)\n');
fprintf(fid, 'the second column shows the number of iterations\n');
fprintf(fid, 'the third column shows the number of symbol errors\n');
fprintf(fid, 'the fourth column shows the number of bit errors\n');
fprintf(fid, 'the sixth column shows SER\n');
fprintf(fid, 'the seventh column shows BER\n');
fprintf(fid, '-------------\n');
fprintf(fid, 'MMSE-OSIC\n');
fprintf(fid, 'SNR \t iteration \t symbol errors \t  bit errors \t SER \t BER\n');



%% Monte-Carlo Simulation
for count=1:length(SNR)    
    symError_sel_MMSE_OSIC=0;  %record the symbol error of MMSE
    bitError_sel_MMSE_OSIC=0; %record the bit error of MMSE
    Realization=0;    %number of channel realization
while(symError_sel_MMSE_OSIC<minSymbolError||Realization<ChannelRealization)
dataIn = randi(M,Nt,1);  % generate the index of transmit bit sequence 
dataMod=zeros(Nt,1); 
grayData=zeros(Nt,1);
for count1=1:Nt      
    dataMod(count1)=symConstell(dataIn(count1));  %generate the modulated symbol vector
    grayData(count1)=graycode(dataIn(count1));   %generate the gray code bit sequence of the modulated symbol vector
end
H=complex(normrnd(0,sqrt(1/2),[Nr,Nt]), normrnd(0,sqrt(1/2),[Nr,Nt]));   %generate channel matrix
n=complex(normrnd(0,sqrt(noiseV(count)/2),Nr,1),normrnd(0,sqrt(noiseV(count)/2),Nr,1));  %generate AWGN vector
sigRec=H*dataMod+n;   %generate receive signal vector
% [symOut_sel_MMSE_OSIC]=sel_MMSE_OSIC(sigRec, H, SNRd(count), M, pav, symConstell);   %Sel_MMSE-OSIC detection
[symOut_sel_MMSE_OSIC]=MMSE_OSIC(sigRec, H, SNRd(count)/Nt, M, pav);   %MMSE-OSIC detection
% [symOut_sel_MMSE_OSIC]=MMSE(sigRec, H, SNRd(count)/Nt, M, pav);   %MMSE detection
symError1_v=abs(dataMod-symOut_sel_MMSE_OSIC);
ErrorIndex=find(symError1_v>1e-5);
symError_sel_MMSE_OSIC=symError_sel_MMSE_OSIC+length(ErrorIndex);  %calculate the symbol error of MMSE in this channel realization
bitOut1=grayDecoder(symOut_sel_MMSE_OSIC, ErrorIndex, graycode, symConstell);  %decoding the gray code of the output symbol vector of MMSE-OSIC
bitError_sel_MMSE_OSIC=bitError_sel_MMSE_OSIC+checkBitError(bitOut1, ErrorIndex, grayData, M);% calculate the bit error of MMSE-OSIC in this channel realization


Realization=Realization+1;  %calculate the number of channel realizations
 end
% Calculate the SER and BER
SER_sel_MMSE_OSIC(count)=symError_sel_MMSE_OSIC/(Realization*Nt);  %caculate symbol error rate of MMSE-OSIC
BER_sel_MMSE_OSIC(count)=bitError_sel_MMSE_OSIC/(Realization*Nt*ceil(log2(M))); %calculate bit error ra%dte of MMSE-OSIC
fprintf(fid,'%d \t %d \t %d \t %d \t %e \t %e\n', SNR(count), Realization, symError_sel_MMSE_OSIC, bitError_sel_MMSE_OSIC, SER_sel_MMSE_OSIC(count), BER_sel_MMSE_OSIC(count));
end
fprintf(fid, '\n');
fprintf(fid, '****************\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'The program ends successfully!\n');
fprintf(fid, '=========================================================================================\n');
fprintf(fid, '\n');
fclose(fid);
toc
