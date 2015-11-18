%% Main System File for Bit Error Rate Monte-Carlo Simulation 
%=============
%Author: Tianpei Chen
% Introduction:
% this routine performs the Monte-Carlo Simulation of the bit error rate of
% proposed detector in Large Scale-MIMO (LS-MIMO) systems. A set of signal
% to noise ratio (SNR) points are considered, for each SNR point, the
% simulations works in the unit of channel realizations, symbol error rate
% and bit error rate of the detectors are calculated when both a certain number
% of symbol errors are accumulated and adequate channel realization times
% is reached.
%Date: Nov 15 2015
%=============
%% Function Description
% grayEncoder: generate gray code list
% symbolConstellation: generate symbol constellation alphabet 
% RSVD: real-Support Vector Detector
% MMSE: minimum mean square error detector 
% grayDecoder: decoding symbol estimations to corresponding gray code
% checkBitError: comparing the output bit sequence and the original
% transmit bit sequence, then calculate bit errors

close all
clear all
tic

%% System Configuration
M =4;         %size of the signal constellation alphabet  (rectangular M-QAM)
Nt=32;         %number of transmit antennas
Nr=32;         %number of receive antennas
SNR=[10:2:26];       %signal to noise ratio (SNR) in (dB)
SNRd=10.^(SNR.*0.1);   %SNR in dicimal
noiseV=1./SNRd;   %noise variance of  Additive White Gaussian Noise (AWGN) 
pav=1/Nt;  %average power of the transmitted symbols
SER_MMSE=zeros(length(SNR),1);    %symbol error rate of MMSE detector for every SNR (dB) 
BER_MMSE=zeros(length(SNR),1);    %bit error rate of MMSE detector for every SNR (dB)
SER_RSVD=zeros(length(SNR),1);    %symbol error rate of RSVD detector for every SNR (dB)
BER_RSVD=zeros(length(SNR),1);    %bit error rate of RSVD detector for every SNR (dB)
outlier_SVD_V=zeros(length(SNR),1);  %number of the extreme outliers of SVR
outlier_MMSE_V=zeros(length(SNR),1); %number of the extreme outliers of MMSE
C=10;              %parameter to control the tradeoff between the penaty term of outliers and regularization term in SVR
epsilon=1e-7;     %parameter to control the estimation precision of SVR
tol=1e-3;         %parameter to control the tolerance of the maximal duality gap (optimization of SVR stops when the duality gap gets smaller than this tolerance)
ChannelRealization=1e3;  %the minimum number of channel realizations
minSymbolError=150;      %the minimum number of symbol errors detected 

%% Signal Modulation
graycode=grayEncoder(M); %gray code encoder
[symConstell]=symbolConstellation( M, pav );  %generate symbol constellation

%% generate file to record the simulation results
fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\Verification Test\Test Data\ExtremelyCase_Test.txt', 'a');
fprintf(fid, '\n');
fprintf(fid, '-----------------\n');
fprintf(fid ,'this file records the BER simulation results of RSVD and MMSE in LS-MIMO systems\n');
fprintf(fid, 'SYSTEM CONFIGURATION\n');
fprintf(fid, '****************\n');
fprintf(fid, '%d X %d MIMO system with %d QAM modulation\n', Nr,Nt,M);
fprintf(fid, 'the average power of transmit symbols is: %0.10f\n', 1/Nt);
fprintf(fid,'the SNR (dB) are:\n');
for count=1:length(SNR)
    fprintf(fid, '%d ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'the hyperparameter settings for RSVD are\n');
fprintf(fid, 'C: %f\n', C);
fprintf(fid, 'tolerance: %e\n', tol);
fprintf(fid, 'epsilon: %e\n', epsilon);
fprintf(fid, 'The minimum channel realization is %e\n', ChannelRealization);
fprintf(fid, 'The minimum symbol error accumulatted is %e\n', minSymbolError);
fprintf(fid, '****************\n');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'Simulation Output\n');
fprintf(fid, '****************\n');
% fprintf(fid, 'BER of RSVD:\n');
fprintf(fid, 'average number of outliers in RSVD\n');


%% Monte-Carlo Simulation
for count=1:length(SNR)    
    symError_RSVD=0;   %record the symbol error of RSVD
    bitError_RSVD=0;   %record the bit error of RSVD
    symError_MMSE=0;  %record the symbol error of MMSE
    bitError_MMSE=0; %record the bit error of MMSE
    Realization=0;    %number of channel realization
    outlier_SVD=0;  %number of the outliers of SVD
    outlier_MMSE=0; %number of the outliers of MMSE
% while(symError_MMSE<minSymbolError||symError_RSVD<minSymbolError||Realization<ChannelRealization)
while(Realization<ChannelRealization)
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
%Transform the complex system model to equivalent real systems model
sigRec_r=[real(sigRec);imag(sigRec)];
H_r=[real(H), -imag(H); imag(H), real(H)];
[symOut_RSVD, outlier_SVD_tmp]=RSVD(H_r,  sigRec_r, SNRd(count),  M, pav, C, tol ,epsilon);  %RSVD detection
[symOut_MMSE, outlier_MMSE_tmp]=MMSE(sigRec, H, SNRd(count), M, pav);   %MMSE detection


outlier_SVD=outlier_SVD+outlier_SVD_tmp;
symError1_v=abs(dataMod-symOut_RSVD);
symError_RSVD=symError_RSVD+length(find(symError1_v>1e-4));  %calculate the symbol error of RSVD in this channel realization
bitOut1=grayDecoder(symOut_RSVD, graycode,symConstell);  %decoding the gray code of the output symbol vector of RSVD
bitError_RSVD=bitError_RSVD+checkBitError(bitOut1, grayData, M);% calculate the bit error of RSVD in this channel realization


outlier_MMSE=outlier_MMSE+outlier_MMSE_tmp;
symError2_v=abs(dataMod-symOut_MMSE);
symError_MMSE=symError_MMSE+length(find(symError2_v>1e-4));  %calculate the symbol error of MMSE in this channel realization
bitOut2=grayDecoder(symOut_MMSE, graycode,symConstell);  %decoding the gray code of the output symbol vector of MMSE
bitError_MMSE=bitError_MMSE+checkBitError(bitOut2, grayData, M);% calculate the bit error of MMSE in this channel realization
Realization=Realization+1;  %calculate the number of channel realizations
 end
% Calculate the SER and BER
SER_RSVD(count)=symError_RSVD/(Realization*Nt);  %caculate symbol error rate of RSVD
BER_RSVD(count)=bitError_RSVD/(Realization*Nt*ceil(log2(M))); %calculate bit error rate of RSVD
SER_MMSE(count)=symError_MMSE/(Realization*Nt);  %caculate symbol error rate of MMSE
BER_MMSE(count)=bitError_MMSE/(Realization*Nt*ceil(log2(M))); %calculate bit error rate of MMSE

% Calculate the average number of the outliers
outlier_SVD_V(count)=outlier_SVD/ChannelRealization;  %calculate the average number of outliers of SVD
outlier_MMSE_V(count)=outlier_MMSE/ChannelRealization; %calculate the average number of outliers of MMSE

% fprintf(fid,'%0.10f, ', BER_RSVD(count));
fprintf(fid, '%0.10f, ', outlier_SVD_V(count) );
end
fprintf(fid, '\n');
% fprintf(fid, 'BER of MMSE:\n');
fprintf(fid, 'average number of outliers in MMSE\n');
for count1=1:length(SNR)
%     fprintf(fid, '%0.10f, ', BER_MMSE(count1));
fprintf(fid, '%0.10f, ', outlier_MMSE_V(count1));
end
fprintf(fid, '\n');
fprintf(fid, '****************\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'The program ends successfully!\n');
fprintf(fid, '-----------------\n');
fclose(fid);
toc