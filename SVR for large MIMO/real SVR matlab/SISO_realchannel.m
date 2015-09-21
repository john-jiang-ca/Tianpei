%SISO channel with M-PAM modulation 
%optimul receiver
M=2; %modulation sheme
Pav=1; %average power per symbol
d=sqrt(3*Pav/((M^2-1)));   %half of the Euclidean distance of two symbols
symbolReal=zeros(M,1);
for i=1:M
    symbolReal(i)=-(M-1)*d+(i-1)*2*d;
end
SNR=[6,8,10];
SNRd=10.^(SNR./10);
NoiseVariance=zeros(1, length(SNR));
NoiseVariance=1./SNRd;
symErrorRate=zeros(1,length(SNR));
for count=1:length(SNR)
    channelRealization=0;
    symError=0;
    while(channelRealization<1e3||symError<200)
data=randi(M,1,1);
CG=normrnd(0, 1, 1,1);   %channel gain 
s=symbolReal(data);
n=normrnd(0, sqrt(NoiseVariance(count)),1,1);
 r=CG*s+n;
% r=s+n;
% r=CG^(-1)*r;
%% ML detection
metric=zeros(1,M);
for count1=1:M
    metric(count1)=norm(r-CG*symbolReal(count1));
end
[value, index]=min(metric);
symOut=symbolReal(index);

 if(abs(s-symOut)>1e-3)
 symError=symError+1;
 end
channelRealization=channelRealization+1;
    end
    symErrorRate(count)=symError/channelRealization;
end





