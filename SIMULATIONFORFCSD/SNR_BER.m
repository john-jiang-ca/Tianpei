%this is for ploting the figure of SNR-BER curve
sim8_4_0_12=[0.126779, 0.055532, 0.015791, 0.002398, 0.000201, 0.000013, 0.000001];
sim8_16_0_16=[0.256291,0.225231,0.189038 ,0.148347,0.094947,0.040334,0.008522,0.000909,0.000060];
sim8_64_0_24=[0.329543, 0.303442, 0.275861, 0.247392, 0.220283, 0.192604, 0.161599, 0.118866, 0.062302, 0.017654, 0.002257, 0.000128, 0.000004];
sim16_4_0_6=[0.032680, 0.005664, 0.000492, 0.000018 ];
sim16_16_0_12=[0.213728, 0.181043, 0.138350, 0.075001, 0.019427, 0.002040, 0.000093];
sim32_4_0_4=[0.002406, 0.000150, 0.000006];
figure (1)
semilogy(0:2:12,sim8_4_0_12,'-*'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
figure (1)
semilogy(0:2:16,sim8_16_0_16,'->'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
figure (1)
semilogy(0:2:24,sim8_64_0_24,'-d'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
figure (1)
semilogy(0:2:6,sim16_4_0_6,'-o'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
figure (1)
semilogy(0:2:12,sim16_16_0_12,'-s'),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
figure (1)
semilogy(0:2:4,sim32_4_0_4,'-p'),xlabel('SNR per bit (dB)'),ylabel('BER'),legend('8X8 4QAM','8X8 16QAM','8X8 64QAM','16X16 4QAM','16X16 16QAM','32X32 4QAM');
hold off;
