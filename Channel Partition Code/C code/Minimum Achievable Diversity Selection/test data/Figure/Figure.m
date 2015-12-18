%% This the m file for the figure plotting 
% for the comparison of maximum diversity channel selection (MDS) and 
% Minimum achievable diversity selection (MAD)
% The performance metric considered here are bit error rate (BER) and
% frame error rate (FER) versus average receive signal to noise ratio (SNR)
% the data come from the corresponding file in the same folder

%% Data 
% BER_MDS_16_4: BER of maximum diversity selection in 16X16 system with 4 QAM
% modulation 
% FER_MAD_16_4_4: FER of minimum achievable diversity selection in 16X16 system with
% 4 QAM modulation, the minimum achievable diversity is 4


BER_MDS_16_4=[0.227928,0.188897,0.141356, 0.0878906 , 0.0364937, 0.00946563, 0.00121562, 8.43769e-05, 2.18343e-06]; %snr 0-16
BER_MAD_16_4_4=[0.228775, 0.189791,0.145853, 0.0975969, 0.0507906, 0.0187031,0.00445312, 0.000572006, 6.73455e-05, 2.90153e-06 ];  %snr 0-18
BER_MAD_16_4_9=[ 0.229984, 0.191062, 0.141838,0.0922656, 0.042625,  0.012875,0.0022, 0.000236072, 9.65151e-06 ];%snr 0-16

BER_MDS_16_16=[0.221552, 0.18933, 0.153348, 0.104134, 0.0416984, 0.00739375,0.000673437, 2.57302e-05 ];%snr 8-22
BER_MAD_16_16_4=[0.223845 , 0.194723, 0.163042, 0.124313,0.0715094, 0.0293125,0.00660625, 0.00141094, 9.24318e-05, 6.72702e-06 ]; %snr 8-26
BER_MAD_16_16_9=[0.222439 , 0.191077, 0.157709, 0.111259, 0.0533094, 0.0146906, 0.00199844, 0.000155307, 5.80553e-06]; %snr 8-24

BER_MAD_32_4_4=[0.229445, 0.191398, 0.149786,0.102313, 0.0562984, 0.0212016, 0.00530156, 0.00070625, 8.49322e-05,6.0447e-06  ]; %snr 0-18
BER_MAD_32_4_9=[0.229477, 0.19105, 0.148106,0.10018, 0.0492078,0.0160859,0.00328281, 0.000368941, 3.70159e-05 ]; %snr 0-16
BER_MAD_32_4_16=[0.229358, 0.190434, 0.146927, 0.0955453, 0.0448656, 0.0129797, 0.00223281, 0.000287571, 1.92706e-05 ]; %snr 0-16





%% Figure plotting
figure (1)
semilogy([0:2:16], BER_MDS_16_4, '-*', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:18], BER_MAD_16_4_4, '-^', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:16], BER_MAD_16_4_9, '-+', 'MarkerSize',5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
title('16X16 system 4QAM');
hold on;
legend('MDS', 'MAD (4)', 'MAD (9)');
hold off;

figure (2)
semilogy([8:2:22], BER_MDS_16_16, '-*', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([8:2:26], BER_MAD_16_16_4, '-^', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([8:2:24], BER_MAD_16_16_9, '-+', 'MarkerSize',5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
title('16X16 system 16QAM');
hold on;
legend('MDS', 'MAD (4)', 'MAD (9)');
hold off;

figure (3)
semilogy([0:2:18], BER_MAD_32_4_4, '-*', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:16], BER_MAD_32_4_9, '-^', 'MarkerSize', 5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:16], BER_MAD_32_4_16, '-+', 'MarkerSize',5), xlabel('SNR (dB)'), ylabel('BER');
hold on;
title('32X32 system 4QAM');
hold on;
legend( 'MAD (4)', 'MAD (9)', 'MAD (16)');
hold off;