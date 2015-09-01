%This is real SVR routine
%author: Tianpei Chen
close all;
clear all;
Nr=100;
Nt=100;
M=4; %modulation sheme
d=sqrt(3/(Nt*(M^2-1)));
symbolReal=zeros(M,1);
for i=1:M
    symbolReal(i)=-(M-1)*d+(i-1)*2*d;
end

SNR=[18];  %signal to noise ratio of output
noiseVariance=zeros(length(SNR)); %noise variance
noiseVariance=1./(10.^(SNR./10));
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

for count=1:length(SNR)
   
    symError=0; %symbol errors
    MSE=0; %average mean square errors
    channel_realization=0;   %channel realization times
    while(symError<10||channel_realization<50)
        %% generate data sample and output (with white gaussian noise)
        transmitNum=randi(M, Nt, 1);
        pH=normrnd(0,1, [Nr, Nt]);
        for j=1:Nt
            x(j)=symbolReal(transmitNum(j));
        end
        n=normrnd(0,sqrt(noiseVariance(count)),Nr, 1);
        y=pH*x+n;      
        %% real SVR 
         
         [symOut, lamida, Theta, G, iteration, MSE]=real_SVR_WSSS1D_2Dsolver(pH, y, 10^(SNR(count)/10), symbolReal, Nr, Nt, M);
%          [symOut, lamida, Theta, G, iteration, MSE]=real_SVR_WSSS1D_2Dsolver_withoutNoise(pH, y, 10^(SNR(count)/10), symbolReal, Nr, Nt, M);
        MSE_original=norm(y-pH*x);
        
        for j=1:Nt
            if(abs(x(j)-symOut(j))>1e-5)
                symError=symError+1;
            end
        end
        channel_realization=channel_realization+1;
    end
    
    SymErrorRate(count)=symError*(1/(channel_realization*Nt));
end


