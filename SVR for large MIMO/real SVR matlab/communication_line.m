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
SymErrorRate=zeros(length(SNR),1);
symOut_candidate=zeros(M^(N),Nt-N);
for count=1:length(SNR)
   
    symError=0; %symbol errors
    MSE=0; %average mean square errors
    channel_realization=0;   %channel realization times
    symError2=0;
    while(symError<50||channel_realization<100)
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
        for count2=1:M^N
            for count3=1:N
                x_candidate(count2,count3)=symbolReal(x_candidate(count2,count3));
            end
        end
        y_candidate=zeros(Nr,M^N);
        symOut_candidate=zeros(Nt,M^N);
        Euclidean=zeros(M^N,1);
        symOut_SVR=zeros(Nt-N,M^N);
        for count2=1:M^(N)
        y_candidate(:,count2)=y-pH(:,1:N)*x_candidate(count2,:)';
        %% real SVR 
         [symOut_SVR(:,count2), lamida, Theta, G, iteration, MSE1]=real_SVR_WSSS1D_2Dsolver(pH(:,N+1:Nt), y_candidate(:,count2), SNRd(count), symbolReal, Nr, Nt-N,...
             M,d);
%          symOut_candidate(count2)=[x_candidate(count2,:),symOut]';
         Euclidean(count2)=MSE1;
         
        end
        [value,index]=min(Euclidean);
        symOut=[x_candidate(index,:),symOut_SVR(:,index)'];
        
%          [symOut, lamida, Theta, G, iteration, MSE]=real_SVR_WSSS1D_2Dsolver_withoutNoise(pH, y, 10^(SNR(count)/10), symbolReal, Nr, Nt, M);
         [  symOut2,   MSE2] = MMSE(pH,  y, SNRd, symbolReal, Nr, Nt, M );
        MSE_original=norm(y-pH*x);
        
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
    
    SymErrorRate(count)=symError*(1/(channel_realization*Nt));
end


