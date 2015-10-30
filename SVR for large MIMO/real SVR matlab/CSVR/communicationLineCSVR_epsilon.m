close all
clear all
tic

% cd('/home/tchen44/Documents/spheredecodingtest/4X4')
%% Signal Modulation and MIMO Channel modeling %% 
M =4;         %size of constellation
Nt=32;         %number of transmit antennas
Nr=32;         %number of receive antennas
% x=6;            %diversity gain that required
SNR=[12];       %signal to noise ratio per bit in dB
SNRd=10.^(SNR.*0.1);   %SNR in dicimal
noiseV=1./SNRd;   %noise variance of AWGN 
BER=zeros(length(SNR),1);         %bit error rate
graycode=grayEncoder(M); %gray encoder
maxGeneration=7;
population=32;
% numFlop=zeros(length(SNR),1);    %complexity of the algorithm
global C;
global tol;
global epsilon;
C=4;
tol =1e-3;
epsilon=1e-7;
SER1=zeros(length(SNR),1);  %symbol error rate of RSVD
BER1=zeros(length(SNR), 1); %bit error rate of RSVD
% SER2=zeros(length(SNR),1); %symbol error rate of MMSE-OSIC
% SER3=zeros(length(SNR),1); %symbol error rate of partial-MIC
% SER4=zeros(length(SNR),1); %symbol error rate of ordered-MIC(ascend)
% SER5=zeros(length(SNR),1); %symbol error rate of ordered-MIC(descend)
% SER6=zeros(length(SNR),1); %symbol error rate of CSVD
% symNum=1e3;   %the number of symbol
% symMap=[11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];    %symbol map
pav=1/Nt;  %average symbol power
[symConstell]=symbolConstellation( M, pav );
% dataInMatrix = reshape(dataIn, [], 2); % Reshape data into binary 4-tuples
% dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers
% dataMod = qammod(dataSymbolsIn, M);
% channelInput = reshape(dataMod, [4, length(dataMod)/4])';    %split data stream into 4 channel
% hMod = comm.RectangularQAMModulator('BitInput', true, ...
%         'ModulationOrder', M, 'NormalizationMethod', 'Average power',...
%  'SymbolMapping','Custom','CustomSymbolMapping',symMap);   %generate a QAM modulation object 
%     hDmod=comm.RectangularQAMDemodulator('BitOutput',true,...
%         'ModulationOrder',M,'NormalizationMethod','Average power',...
%     'SymbolMapping','Custom','CustomSymbolMapping',symMap); %generate a QAM demodulation object

%      channelInput = reshape(dataMod, [],Nt);    %split data stream into 4 channel
% hMIMO = comm.MIMOChannel('TransmitCorrelationMatrix',eye(Nt),...
%     'ReceiveCorrelationMatrix',eye(Nr),...
%              'RandomStream', 'mt19937ar with seed',...
%        'PathGainsOutputPort', true);    %MIMO channel modeling
%    dataMod1=dataMod*sqrt(10);
%    dataMod2=[real(dataMod1);imag(dataMod1)];
% hBER = comm.ErrorRate;
% ber = step(hBER, dataIn, double(dataMod2(:)));
% [~, pathGain] = step(hMIMO, channelInput);  %get the receive vectors through MIMO channel
% BitTable = de2bi(symMap, log2(M), 'left-msb');
% hSpDec = comm.SphereDecoder('Constellation', constellation(hMod),...
%         'BitTable', BitTable, 'DecisionType', 'Hard');
%     hBER = comm.ErrorRate;
 fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\test data\BER_RSVD_GA.txt', 'a');
 fprintf(fid, '\n');
 fprintf(fid, '-----------------\n');
 fprintf(fid ,'this file record the simulation results of RSVD-GA\n');
 fprintf(fid, 'the population is %d\n', population);
 fprintf(fid, 'the number of generation is %d\n', maxGeneration);
 fprintf(fid, 'RSVD-GA with OMIC\n');
 fprintf(fid, 'the order method 1\n');
fprintf(fid, 'this is the bit error rate of %d X %d MIMO system with %d QAM modulation\n', Nr,Nt,M);
fprintf(fid, 'the hyperparameters for RSVD-OMIC are\n');
fprintf(fid, 'test the epsilon and C\n');
% fprintf(fid, 'epsilon: %0.10f\n', epsilon);
fprintf(fid, 'tolerence: %0.10f\n', tol);
fprintf(fid, 'without the clipper\n');
% fprintf(fid, 'use the abs of feasibility gap\n');
fprintf(fid,'the SNR are:\n');
for count=1:length(SNR)
    fprintf(fid, '%d ', SNR(count));
end
fprintf(fid, '\n');

fprintf(fid,'the C are:\n');
for count=1:length(C)
    fprintf(fid, '%d ', C(count));
end
fprintf(fid, '\n');

fprintf(fid,'the epsilon are:\n');
for count=1:length(epsilon)
    fprintf(fid, '%d ', epsilon(count));
end
fprintf(fid, '\n');
fprintf(fid,'the symbol error rates are\n');
fclose(fid);
%% Monte-Carlo simulation
for count=1:length(SNR)     %under the SNR from 0 to 10
%     symError1=zeros(length(epsilon),length(C));
    symError1=0;
    symError2=0;
    symError3=0;
    symError4=0;
    symError5=0;
    symError6=0;
    bitError1=0;
    channelRealization=0;          %number of channel realization
%     bitOutput=zeros(nBits,1);
 fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\test data\BER_RSVD_GA.txt', 'a');
 while(symError1<100||channelRealization<1e3)
% symOut=zeros(Nt,1);
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
    grayData(count1)=graycode(dataIn(count1));
end
H=complex(normrnd(0,sqrt(1/2),[Nr,Nt]), normrnd(0,sqrt(1/2),[Nr,Nt]));   %channel matrix
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
  maxStage1=4;
  maxStage2=6;
  symOut1=zeros(Nt,1);
  W=ones(Nt,Nt);
  G=zeros(Nt,Nr);
%   [symOut1]=MIC_Recursive(sigRec, H, symOut1, SNRd(count), M, pav, stage,maxStage, W,G, symConstell);
  symOut_total=zeros(Nt,1);
%   tol=0.1;
  list=1:Nt;
  W=zeros(Nt,1);
  sym_prev=zeros(Nt,1);
  sym_total=zeros(Nt,1);
  G=zeros(Nt);
  order=zeros(Nt,1);
%   [ symOut1 ] = MIC_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell);
%   [symOut2]=MMSE_OSIC(sigRec, H, SNRd(count), M, pav);
%    [symOut3]=Partial_MIC_recursive( sigRec, H, list, W, sym_prev, sym_total, SNRd(count), M, pav, stage, maxStage2, symConstell, tol);
% [ symOut4 ] = MIC_ordering_ascend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order);
% [ symOut1 ] = MIC_ordering_descend_Recursive(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order);
% [symOut5]=CSVD(sigRec, H, SNRd(count), M);
sigRec_r=[real(sigRec);imag(sigRec)];
H_r=[real(H), -imag(H); imag(H), real(H)];
% symOut=zeros(Nt,length(C));
% for count0=1:length(epsilon)
% for count1=1:length(C)
%  [symOut1]=RSVR(H_r,  sigRec_r, SNRd(count),  M, pav, C, tol ,epsilon);
% [ symOut1 ] = RSVR_OMIC_MMSE(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order,C, tol, epsilon);
% [ symOut1 ] = RSVR_OMIC_total(sigRec, H, sym_prev, SNRd(count), M, pav, stage, maxStage1,W,G, symConstell, order,C, tol, epsilon);
generation=0;

symGeneration=zeros(Nt, population);
[ symOut1 ] = RSVD_GA( sigRec, H, symGeneration, M, symConstell, C, tol, epsilon, SNRd(count), pav, population, generation, maxGeneration );
 symError_v=abs(dataMod-symOut1);
symError1=symError1+length(find(symError_v>1e-4));
bitOut1=grayDecoder(symOut1, graycode,symConstell);
bitError1=bitError1+checkBitError(bitOut1, grayData, M);
% end
% end
%  bitOut=step(hDmod,symOut);
%      rxBits = step(hSpDec, channelOutput, squeeze(H));
%      ber = step(hBER, dataIn, double(rxBits(:)));
% for j=1:Nt
% %     if(abs(dataMod(j)-symOut1(j))>1e-4)
% %         symError1=symError1+1;
% %     end
% %     if(abs(dataMod(j)-symOut2(j))>1e-4)
% %         symError2=symError2+1;
% %     end
% %         if(abs(dataMod(j)-symOut3(j))>1e-4)
% %         symError3=symError3+1;
% %         end
% %         if(abs(dataMod(j)-symOut4(j))>1e-4)
% %         symError4=symError4+1;
% %         end
%          if(abs(dataMod(j)-symOut5(j))>1e-4)
%         symError5=symError5+1;
%          end
%           if(abs(dataMod(j)-symOut6(j))>1e-4)
%         symError6=symError6+1;
%         end
% end

% for j=1:nBits
%     if(abs(dataIn(j)-bitOut(j))>1e-4)
%         bitError=bitError+1;
%     end
% end

channelRealization=channelRealization+1;
 end
%  for count0=1:length(epsilon)
%  for count1=1:length(C)
SER1(count)=symError1/(channelRealization*Nt);  %caculate symbol error rate of MIC
BER1(count)=bitError1/(channelRealization*Nt*ceil(log2(M))); %calculate bit error rate of RSVD
fprintf(fid,'%0.10f, ', BER1(count));
%  end
%  fprintf(fid, '\n');
% end
% fprintf(fid, '\n');
% SER2(count)=symError2/(channelRealization*Nt);  %calculate symbol error rate of MMSE-OSIC
%  SER3(count)=symError3/(channelRealization*Nt); %calculate symbol error rate of partial-MIC
% SER4(count)=symError4/(channelRealization*Nt); %calculate symbol error rate of ordered-MIC(ascend)
% SER5(count)=symError5/(channelRealization*Nt); %calculate symbol error rate of ordered-MIC(descend)
% SER6(count)=symError6/(channelRealization*Nt); %calculate symbol error rate of Complex Support Vector Detector(CSVD)
% BER(count)=bitError/(channelRealization*nBits);   %caculate bit error rate
end
fprintf(fid,'\n');
fprintf(fid, 'The program ends successfully!\n');
fprintf(fid, '-----------------\n');
fclose(fid);
% fid1=fopen('BER_8.txt','w');
% fprintf(fid1,'%e\n',BER);
% fclose(fid1);
% fid2=fopen('SER_8.txt','w');
% fprintf(fid2,'%e\n',symErrorrate);
% fclose(fid2);

toc