close all
clear all
% cd('/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/MMSE-OSIC/32/4QAM')
tic
% cd('/home/tchen44/Documents/spheredecodingtest/4X4')
%% Signal Modulation and MIMO Channel modeling %% 
M =4;         %size of constellation
Nt=64;         %number of transmit antennas
Nr=64;         %number of receive antennas
% x=6;            %diversity gain that required
SNR=[10:2:24];       %signal to noise ratio per bit in dB
SNRd=10.^(SNR.*0.1);   %SNR in dicimal
noiseV=1./SNRd;   %noise variance of AWGN 
BER=zeros(length(SNR),1);         %bit error rate
% numFlop=zeros(length(SNR),1);    %complexity of the algorithm
SER1=zeros(length(SNR),1);  %symbol error rate of MMSE-OSIC
% SER2=zeros(length(SNR),1); %symbol error rate of MMSE-OSIC
% SER3=zeros(length(SNR),1); %symbol error rate of partial-MIC
% SER4=zeros(length(SNR),1); %symbol error rate of ordered-MIC(ascend)
% SER5=zeros(length(SNR),1); %symbol error rate of ordered-MIC(descend)
% SER6=zeros(length(SNR),1); %symbol error rate of CSVD
% SER7=zeros(length(SNR),1); %symbol error rate of RSVD aided ordered MIC
pav=1/Nt;  %average symbol power
[symConstell]=symbolConstellation( M, pav );

%% generate files
fid_t=fopen('/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/MMSE-OSIC/64/4QAM/data/README.txt','a');   %the file for explanations
fprintf(fid_t, 'this is %d X %d MIMO\n', Nr, Nt);
fprintf(fid_t, '%d QAM modulation\n', M);
fprintf('\n');
fprintf(fid_t, 'The detectors tested are\n');
fprintf(fid_t, 'MMSE-OSIC-SER1\n');
% fprintf(fid_t, 'RSVD-SER2\n');
% fprintf(fid_t, 'MMSE-OSIC-SER3\n');
% fprintf(fid_t, 'MMSE-OMIC(ascend)-SER4\n');
% fprintf(fid_t, 'MMSE-OMIC(descend)-SER5\n');
% fprintf(fid_t, 'RSVD-OMIC(descend)-SER6\n');


fid1=fopen('/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/MMSE-OSIC/64/4QAM/data/SER1.txt','a');  %the file store the SER values
% fid2=fopen('SER2.txt','a');
% fid3=fopen('SER3.txt','a');
% fid4=fopen('SER4.txt','a');
% fid5=fopen('SER5.txt','a');
% fid6=fopen('SER6.txt','a');
% fid7=fopen('SER6.txt','a');
fprintf(fid1, 'symbol error rate of MMSE-OSIC\n');
% fprintf(fid2, 'symbol error rate of RSVD\n');
% fprintf(fid3, 'symbol error rate of MMSE-OSIC\n');
% fprintf(fid4, 'symbol error rate of MMSE-OMIC(ascend)\n');
% fprintf(fid5, 'symbol error rate of MMSE-OMIC(descend)\n');
% fprintf(fid6, 'symbol error rete of RSVD-OMIC(descend)\n');

fprintf(fid_t, 'SNR in dB are\n');
fprintf(fid1, 'SNR in dB are\n');
% fprintf(fid2, 'SNR in dB are\n');
% fprintf(fid3, 'SNR in dB are\n');
% fprintf(fid4, 'SNR in dB are\n');
% fprintf(fid5, 'SNR in dB are\n');
% fprintf(fid6, 'SNR in dB are\n');
% fprintf(fid7, 'SNR in dB are\n');
for count=1:length(SNR)
    fprintf(fid_t, '%d ', SNR(count));
    fprintf(fid1, '%d ', SNR(count));
%     fprintf(fid2, '%d ', SNR(count));
%     fprintf(fid3, '%d ', SNR(count));
%     fprintf(fid4, '%d ', SNR(count));
%     fprintf(fid5, '%d ', SNR(count));
%     fprintf(fid6, '%d ', SNR(count));
%     fprintf(fid7, '%d ', SNR(count));
end
fprintf(fid_t, '\n');
fprintf(fid1, '\n');
% fprintf(fid2, '\n');
% fprintf(fid3, '\n');
% fprintf(fid4, '\n');
% fprintf(fid5, '\n');
% fprintf(fid6, '\n');
% fprintf(fid7, '\n');

fprintf(fid_t, 'SNR in decimal are\n');
fprintf(fid1, 'SNR in decimal are\n');
% fprintf(fid2, 'SNR in decimal are\n');
% fprintf(fid3, 'SNR in decimal are\n');
% fprintf(fid4, 'SNR in decimal are\n');
% fprintf(fid5, 'SNR in decimal are\n');
% fprintf(fid6, 'SNR in decimal are\n');
% fprintf(fid7, 'SNR in decimal are\n');
for count=1:length(SNRd)
        fprintf(fid_t, '%f ', SNRd(count));
    fprintf(fid1, '%f ', SNRd(count));
%     fprintf(fid2, '%f ', SNRd(count));
%     fprintf(fid3, '%f ', SNRd(count));
%     fprintf(fid4, '%f ', SNRd(count));
%     fprintf(fid5, '%f ', SNRd(count));
%     fprintf(fid6, '%f ', SNRd(count));
%     fprintf(fid7, '%f ', SNRd(count));
end

fprintf(fid_t, '\n');
fprintf(fid1, '\n');
% fprintf(fid2, '\n');
% fprintf(fid3, '\n');
% fprintf(fid4, '\n');
% fprintf(fid5, '\n');
% fprintf(fid6, '\n');
% fprintf(fid7, '\n');

fclose(fid_t);
fclose(fid1);
% fclose(fid2);
% fclose(fid3);
% fclose(fid4);
% fclose(fid5);
% fclose(fid6);
% fclose(fid7);



 
%% Monte-Carlo simulation
for count=1:length(SNR)     %under the SNR from 0 to 10
    symError1=0;
    symError2=0;
    symError3=0;
    symError4=0;
    symError5=0;
    symError6=0;
    symError7=0;
    bitError=0;
    channelRealization=0;          %number of channel realization
%     bitOutput=zeros(nBits,1);

fid1=fopen('/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/MMSE-OSIC/64/4QAM/data/SER1.txt','a');  %the file store the SER values
% fid2=fopen('SER2.txt','a');
% fid3=fopen('SER3.txt','a');
% fid4=fopen('SER4.txt','a');
% fid5=fopen('SER5.txt','a');
% fid6=fopen('SER6.txt','a');
% fid7=fopen('SER6.txt','a');
 while(symError1<200||channelRealization<5*1e3)
% symOut=zeros(Nt,1)
% for k=1:symNum/Nt
% for i=1:Nt
%     for j=1:Nr
%         A(i,j)=pathGain(k,1,i,j);
%     end
% end
% H(k,1,:,:)=A;
% H=A';
% nBits =Nt*log2(M);   %;number of bits transmit
dataIn = randi(M,Nt,1);  % Generate vector of input data (1 to M)
% dataMod = step(hMod, dataIn);  %modulate the bit sequence
dataMod=zeros(Nt,1);
for count1=1:Nt
    dataMod(count1)=symConstell(dataIn(count1));
end
H=complex(normrnd(0,1/2,[Nr,Nt]), normrnd(0,1/2,[Nr,Nt]));   %channel matrix
%      H(k,1,:,:)=H1';
% channelOutput=H*dataMod;
% hAWGN = comm.AWGNChannel('EbNo',SNR(i),'SignalPower',10,...
%         'BitsPerSymbol',log2(M));
% hAWGN=comm.AWGNChannel('NoiseMethod','Variance',...
%     'VarianceSource','Property','Variance',noiseVariance(i)/2);
n=complex(normrnd(0,sqrt(noiseV(count)/2),Nr,1),normrnd(0,sqrt(noiseV(count)/2),Nr,1));
%     ,'VarianceSource', 'Property', 'Variance', noiseVariance(i));   %generate AWGN noise object
%     sigReceive(k,1:Nr)=step(hAWGN, channelOutput);   %add AWGN to channelOutput
%     [ H, P] = Ordering(H);   %apply channel ordering
% sigReceive(k,:)=channelOutput;
sigRec=H*dataMod+n;
%      sigVector=sigReceive(k,:)';
%   [ symOutput(Nt*(k-1)+1:Nt*k),N,d ] = sel_MMSE( sigVector,H,x,(10.^(SNR.*0.1)),dataMod);    %detection using sel-MMSE algorithm
%[ symOutput1(Nt*(k-1)+1:Nt*k,1),numV,numF,numL] = SphDec_chen(H,sigReceive(k,:),dataMod(1:Nt,1));
%  symOutput(Nt*(k-1)+1:Nt*k,1)=P*symOutput(Nt*(k-1)+1:Nt*k,1);    %change the order of transmit symbol vector back

% numFlop(i)=numFlop(i)+numF;
% end
% symOut=dataMod;

  stage=1;
%   maxStage0=2;
  maxStage1=4;
%   maxStage2=8;
  symOut1=zeros(Nt,1);
  W=ones(Nt,Nt);
  G=zeros(Nt,Nr);
%   [symOut1]=MIC_Recursive(sigRec, H, symOut1, SNRd(count), M, pav, stage,maxStage, W,G, symConstell);
  symOut_total=zeros(Nt,1);
  tol=0.1;
  list=1:Nt;
  W=zeros(Nt,1);
  sym_prev=zeros(Nt,1);
  sym_total=zeros(Nt,1);
  G=zeros(Nt);
  order=zeros(Nt,1);
  H_r=[real(H), -imag(H); imag(H), real(H)];
sigRec_r=[real(sigRec);imag(sigRec)];
%   [symOut1]=MMSE_complex(sigRec, H, SNRd(count), M, pav);
%   [ symOut1 ] = MIC_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell);
%    [symOut2]=RSVR(H_r,  sigRec_r, SNRd(count),  M, pav);
  [symOut1]=MMSE_OSIC(sigRec, H, SNRd(count), M, pav);

%    [symOut3]=Partial_MIC_recursive( sigRec, H, list, W, sym_prev, sym_total, SNRd(count), M, pav, stage, maxStage2, symConstell, tol);
% [symOut3]=RSVD_OSIC(sigRec, H, SNRd(count), M, pav);
% [ symOut4 ] = MIC_ordering_ascend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order);
% [ symOut4 ] = MIC_ordering_descend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage0,W,G, symConstell, order);
% [ symOut5 ] = MIC_ordering_descend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order);
% [ symOut6 ] = MIC_ordering_descend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage2,W,G, symConstell, order);
% [symOut5]=CSVD(sigRec, H, SNRd(count), M);
% [ symOut6 ] = RSVR_MIC_ordering_descend(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order);

%  bitOut=step(hDmod,symOut);
%      rxBits = step(hSpDec, channelOutput, squeeze(H));
%      ber = step(hBER, dataIn, double(rxBits(:)));
% for j=1:Nt
%     if(abs(dataMod(j)-symOut1(j))>1e-4)
%         symError1=symError1+1;
%     end
%     if(abs(dataMod(j)-symOut2(j))>1e-4)
%         symError2=symError2+1;
%     end
%         if(abs(dataMod(j)-symOut3(j))>1e-4)
%         symError3=symError3+1;
%         end
%         if(abs(dataMod(j)-symOut4(j))>1e-4)
%         symError4=symError4+1;
%         end
%          if(abs(dataMod(j)-symOut5(j))>1e-4)
%         symError5=symError5+1;
%          end
%           if(abs(dataMod(j)-symOut6(j))>1e-4)
%         symError6=symError6+1;
%           end
%        if(abs(dataMod(j)-symOut7(j))>1e-4)
%         symError7=symError7+1;
%         end
% end
symError1_v=abs(dataMod-symOut1);
symError1_v=find(symError1_v>1e-4);
symError1=symError1+length(symError1_v);
% for j=1:nBits
%     if(abs(dataIn(j)-bitOut(j))>1e-4)
%         bitError=bitError+1;
%     end
% end

channelRealization=channelRealization+1;
 end
SER1(count)=symError1/(channelRealization*Nt);  %caculate symbol error rate of MMSE-OSIC
% SER2(count)=symError2/(channelRealization*Nt);  %calculate symbol error rate of RSVD
% SER3(count)=symError3/(channelRealization*Nt); %calculate symbol error rate of MMSE-OSIC
% SER4(count)=symError4/(channelRealization*Nt); %calculate symbol error rate of ordered-MIC(ascend)
% SER5(count)=symError5/(channelRealization*Nt); %calculate symbol error rate of ordered-MIC(descend)
% SER6(count)=symError6/(channelRealization*Nt); %calculate symbol error rate of RSVD-OMIC(descend)
fprintf(fid1, '%0.10f ', SER1(count));
% fprintf(fid2, '%0.10f ', SER2(count));
% fprintf(fid3, '%0.10f ', SER3(count));
% fprintf(fid4, '%0.10f ', SER4(count));
% fprintf(fid5, '%0.10f ', SER5(count));
% fprintf(fid6, '%0.10f ', SER6(count));
% fprintf(fid7, '%0.10f ', SER7(count));
fclose(fid1);
% fclose(fid2);
% fclose(fid3);
% fclose(fid4);
% fclose(fid5);
% fclose(fid6);
% fclose(fid7);
% SER7(count)=symError7/(channelRealization*Nt); %calculate symbol error rate of Complex Support Vector Detector(CSVD)
% BER(count)=bitError/(channelRealization*nBits);   %caculate bit error rate
end
fid_t=fopen('/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/MMSE-OSIC/64/4QAM/data/README.txt', 'a');
fprintf(fid_t, 'The whole program end successfully!');
fclose(fid_t);
% fid1=fopen('BER_8.txt','w');
% fprintf(fid1,'%e\n',BER);
% fclose(fid1);
% fid2=fopen('SER_8.txt','w');
% fprintf(fid2,'%e\n',symErrorrate);
% fclose(fid2);

toc
