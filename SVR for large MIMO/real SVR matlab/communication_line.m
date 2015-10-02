%This is real SVR routine
%author: Tianpei Chen
close all;
clear all;
cd 'F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab'
Nr=32;
Nt=32;
N=4;
M=4; %modulation sheme
d=sqrt(3/(Nt*(M^2-1)));
symbolReal=zeros(M,1);
epsilon=2*(1e-1);
C=1;
tol=1e-2;
tau=0.2;
for i=1:M
    symbolReal(i)=-(M-1)*d+(i-1)*2*d;
end

SNR=[20];  %signal to noise ratio of output
noiseVariance=zeros(length(SNR)); %noise variance
SNRd=zeros(length(SNR),1);
SNRd=10.^(SNR./10);
noiseVariance=1./(SNRd);
pH=zeros(Nr,Nt); %data samples 
y=zeros(Nr,1);   %output
x=zeros(Nt,1);   %real coefficience  
n=zeros(Nr,1); %
lamida=zeros(Nr); %dual variable pair
K=zeros(Nr, Nr); %kernel matrix
Phi=zeros(Nr); %intermediate variables
iteration=0; %iteration time of SVR algorithm
symError=0;
SymErrorRate1=zeros(length(SNR),1);
SymErrorRate2=zeros(length(SNR),1);
symOut_candidate=zeros(M^(N),Nt-N);
%% document settings
fid1=fopen('D_SVR.txt','a');
fid2=fopen('D_MMSE.txt','a');
fid3=fopen('D_original.txt','a');

fprintf(fid1, 'This document stores the information for SVR\n');
fprintf(fid1, '%d X %d MIMO system\n',Nr,Nt);
fprintf(fid1, '%d PAM modulation scheme\n', M);
fprintf(fid1, 'SNR are\n');
fprintf(fid1, '   ');
for i=1:length(SNR)
    fprintf(fid1, '%d  ', SNR(i));
end
fprintf(fid1, '\n\n');
fprintf(fid1, '   ');
for i=1:length(SNR)
    fprintf(fid1, '%f  ', SNRd(i));
end
fprintf(fid1, '\n\n');
fprintf(fid1, 'SER AMMSE\n');
fprintf(fid1, 'hyperparameters\n');
fprintf(fid1, 'epsilon is %0.10f\n',epsilon);
fprintf(fid1,'C is %f\n', C);
fprintf(fid1, 'tolerance is %0.8f\n', tol);
fclose(fid1);


fprintf(fid2, 'This document stores the information for MMSE\n');
fprintf(fid2, '%d X %d MIMO system\n',Nr,Nt);
fprintf(fid2, '%d PAM modulation scheme\n', M);
fprintf(fid2, 'SNR are\n');
fprintf(fid2, '   ');
for i=1:length(SNR)
    fprintf(fid2, '%d  ', SNR(i));
end
fprintf(fid2, '\n\n');
fprintf(fid2, '   ');
for i=1:length(SNR)
    fprintf(fid2, '%f  ', SNRd(i));
end
fprintf(fid2, '\n\n');
fprintf(fid2, 'SER AMMSE\n');
fprintf(fid2, 'hyperparameters\n');
fprintf(fid2, 'epsilon is %0.10f\n',epsilon);
fprintf(fid2,'C is %f\n', C);
fprintf(fid2, 'tolerance is %0.8%f\n', tol);
fclose(fid2);


fprintf(fid3, 'This document stores the information for original MMSE\n');
fprintf(fid3, '%d X %d MIMO system\n',Nr,Nt);
fprintf(fid3, '%d PAM modulation scheme\n', M);
fprintf(fid3, 'SNR are\n');
fprintf(fid3, '  ');
for i=1:length(SNR)
    fprintf(fid3, '%d  ', SNR(i));
end
fprintf(fid3, '\n\n');
fprintf(fid3, '  ');
for i=1:length(SNR)
    fprintf(fid3, '%f  ', SNRd(i));
end
fprintf(fid3, '\n\n');
fprintf(fid3, 'AMMSE\n');
fprintf(fid3, 'hyperparameters\n');
fprintf(fid3, 'epsilon is %0.10f\n',epsilon);
fprintf(fid3,'C is %f\n', C);
fprintf(fid3, 'tolerance is %0.8%f\n', tol);
fclose(fid3);

%% channel realizations
for count=1:length(SNR)
   fid1=fopen('D_SVR.txt','a');
fid2=fopen('D_MMSE.txt','a');
fid3=fopen('D_original.txt','a');
    symError=0; %symbol errors
    MSE_SUM1=0; %average mean square error of SVR
    MSE_SUM2=0; %average mean square error of MMSE
    MMSE_SUM=0; %average minimum mean square error of original one 
    channel_realization=0;   %channel realization times
    symError2=0;
    while(symError<300||channel_realization<1e2)
        %% generate data sample and output (with white gaussian noise)
        transmitNum=randi(M, Nt, 1);
        pH=normrnd(0,1, [Nr, Nt]);
        for j=1:Nt
            x(j)=symbolReal(transmitNum(j));
        end
        n=normrnd(0,sqrt(noiseVariance(count)),Nr, 1);
        y=pH*x+n;    
        
        
        %% channel partition
     
         x_candidate=zeros(M^(N),N);
        x_candidate=fullfact([M M M M]);
        x_candidate=x_candidate';
        for count2=1:M^N
            for count3=1:N
                x_candidate(count3,count2)=symbolReal(x_candidate(count3,count2));
            end
        end
        y_candidate=zeros(Nr,M^N);
        symOut_candidate=zeros(Nt,M^N);
        Euclidean=zeros(M^N,1);
        symOut_SVR=zeros(Nt-N,M^N);
%         for count2=1:M^(N)
%          y_candidate(:,count2)=y-pH(:,1:N)*x_candidate(:,count2);
        %% real SVR 
%           [symOut_SVR(:,count2), lamida, Theta, G, iteration, MSE1]=real_SVR_WSSS1D_2Dsolver(pH(:,N+1:Nt), y_candidate(:,count2), SNRd(count), symbolReal, Nr, Nt-N,...
%               M,d);
[symOut, lamida, Theta, G, iteration, MSE1, NumReliable]=real_SVR_WSSS1D_2Dsolver(pH, y, SNRd(count), symbolReal, Nr, Nt, M,d,epsilon, C, tol,tau);
numberR=0;
subsymOut1=[];
subpH1=[];
subpH2=[];
for j=1:Nt
  numberR=numberR+NumReliable(j);
  if(NumReliable(j)==1)
      subsymOut1=[subsymOut1;symOut(j)];
      subpH1=[subpH1,pH(:,j)];
  else
      subpH2=[subpH2,pH(:,j)];
  end
end
y1=zeros(Nr,1);
y1=y-subpH1*subsymOut1;
subsymOut2=zeros(Nt-numberR,1);
[subsymOut2, lamida, Theta, G, iteration, MSE1, NumReliable2]=real_SVR_WSSS1D_2Dsolver(subpH2, y1, SNRd(count), symbolReal, Nr, Nt-numberR, M,d,epsilon, C, tol,tau);
index1=0;
index0=0;
for j=1:Nt
    if(NumReliable(j)==1)
        index1=index1+1;
        symOut(j)=subsymOut1(index1);
    else
        index0=index0+1;
        symOut(j)=subsymOut2(index0);
    end
end



    
%           symOut_candidate(count2)=[x_candidate(:,count2);symOut'];
%           Euclidean(count2)=MSE1;
         
%         end
%         [value,index]=min(Euclidean);
%         symOut=[x_candidate(index,:),symOut_SVR(:,index)'];
        
%          [symOut, lamida, Theta, G, iteration, MSE]=real_SVR_WSSS1D_2Dsolver_withoutNoise(pH, y, 10^(SNR(count)/10), symbolReal, Nr, Nt, M);
        MSE_SUM1=MSE_SUM1+MSE1;
         [  symOut2,   MSE2] = MMSE(pH,  y, SNRd(count), symbolReal, Nr, Nt, M );
        MSE_SUM2=MSE_SUM2+MSE2;
% symOut2=zeros(Nt,1);
        MSE_original=norm(y-pH*x);
        MMSE_SUM=MMSE_SUM+MSE_original;
        
        for j=1:Nt
            if(abs(x(j)-symOut(j))>1e-5)
                symError=symError+1;
            end
        end
        for m=1:Nt
            if(abs(x(m)-symOut2(m))>1e-5)
                symError2=symError2+1;
            end
        end
        channel_realization=channel_realization+1;
    end
    MSE_SUM1=MSE_SUM1/(channel_realization);
    MSE_SUM2=MSE_SUM2/(channel_realization);
    MMSE_SUM=MMSE_SUM/(channel_realization);
   SymErrorRate1(count)=symError*(1/(channel_realization*Nt));
   SymErrorRate2(count)=symError2*(1/(channel_realization*Nt));
   fprintf(fid1, 'symbol error rate is %0.10f ', SymErrorRate1(count));
   fprintf(fid1, 'Average Mean Square Error is %0.10f ', MSE_SUM1);
   fprintf(fid2, 'symbol error rate is %0.10f ', SymErrorRate2(count));
   fprintf(fid2, 'Average Mean Square Error is %0.10f ', MSE_SUM2);
   fprintf(fid3,'original Mean Square Error is %0.10f ', MMSE_SUM);
   fclose(fid1);
   fclose(fid2);
   fclose(fid3);
end


