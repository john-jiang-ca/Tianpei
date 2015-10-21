%this is for ploting the figure of SNR-BER curve
sim8_4_0_12=[0.126779, 0.055532, 0.015791, 0.002398, 0.000201, 0.000013, 0.000001];
sim8_16_0_16=[0.256291,0.225231,0.189038 ,0.148347,0.094947,0.040334,0.008522,0.000909,0.000060];
sim8_64_0_24=[0.329543, 0.303442, 0.275861, 0.247392, 0.220283, 0.192604, 0.161599, 0.118866, 0.062302, 0.017654, 0.002257, 0.000128, 0.000004];
sim16_4_0_6=[0.032680, 0.005664, 0.000492, 0.000018 ];
sim16_16_0_12=[0.213728, 0.181043, 0.138350, 0.075001, 0.019427, 0.002040, 0.000093];
0.21539718032 0.13699733913 0.06087859422 0.02043750092,0.002406, 0.000150, 0.000006
sim32_4_0_4=[0.002406, 0.000150, 0.000006];
sim48_4_0_4=[0.000511, 0.000032, 0.000000548];
0.25442603230, 0.19931979477, 0.12637500465,0.05224, 0.01231614593 ,0.000511, 0.000032, 0.000000548
fig1=figure;
Xsize=8.5;Ysize=11;
% Xleft=(8-Xsize)/2;Ytop=(12-Ysize)/2;
Xleft=0;
Ytop=0;
set(fig1,'Units','inches');
set(fig1,'position', [Xleft,Ytop,Xsize,Ysize]);
set(fig1,'PaperUnits','inches');
set(fig1,'PaperPosition',[Xleft,Ytop,Xsize,Ysize]);
set(gca,'fontsize',11);
fig1=semilogy(0:2:12,sim8_4_0_12,'-*');xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:16,sim8_16_0_16,'->'),xlabel('SNR per bit (dB)'),ylabel('BER');
fig1=semilogy(0:2:24,sim8_64_0_24,'-d'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:6,sim16_4_0_6,'-o'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:12,sim16_16_0_12,'-s'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:4,sim32_4_0_4,'-p'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:4,sim48_4_0_4,'-+'),xlabel('SNR per bit (dB)'),ylabel('BER'),legend('8X8 4QAM','8X8 16QAM','8X8 64QAM','16X16 4QAM','16X16 16QAM','32X32 4QAM','48X48 4QAM');
hold on;
print('-depsc','-tiff','-r300','BER_curves');
hold off;

