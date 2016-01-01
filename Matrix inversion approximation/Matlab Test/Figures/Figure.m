%% This the m file for the figure plotting 
% The test of the performances of linear detectors (LD) (MMSE) that exploiting Neumann series expansion 1-Newton Iteration
% Approximation (SENIA-1)
% MMSE with exact matrix inversion (EMI) is also considered as comparison
% The performance metric considered here are bit error rate (BER) versus average receive signal to noise ratio (SNR)
% the data come from the corresponding file in the same folder

%% Data example
% BER_E_128_8_16: BER of MMSE exploiting exact matrix inversion for 128X8
% system with 16QAM modulation
% BER_A_128_8_16: BER of MMSE exploiting SENIA-1 matrix inversion

%% 16QAM
%The expectation of OM is 0.79999356
BER_E_128_8_16=[3.365625e-02, 1.137500e-02,2.956481e-03,2.557492e-04,8.469799e-06]; %snr 0-8
BER_A_128_8_16=[3.456250e-02,1.153125e-02,2.986045e-03, 2.915541e-04,1.058725e-05 ];  %snr 0-8

%The expectation of OM is 1.442349e-02
BER_E_128_32_16=[1.738750e-01,1.280547e-01,8.495313e-02,4.814844e-02,1.981250e-02,6.234375e-03,...
    8.750000e-04,4.112059e-05,7.034182e-07]; %snr 0-16
BER_A_128_32_16=[1.855859e-01, 1.506328e-01,1.253672e-01,1.051875e-01,8.781250e-02,7.839844e-02,6.975781e-02,...
    6.941813e-02,6.686763e-02];  %snr 0-16

%The expectation of OM is 0.63402790
BER_E_64_8_16=[9.053125e-02,5.409375e-02,2.393750e-02,7.093750e-03,1.249011e-03,1.131467e-04,2.449150e-06]; %snr 0-12
BER_A_64_8_16=[9.365625e-02,5.718750e-02,2.718750e-02,8.906250e-03,2.473288e-03,5.510247e-04,1.884621e-04];  %snr 0-12

%The expectation of OM is 1.800414e-13
BER_E_32_32_16=[8.879688e-03,5.419531e-03,3.740625e-03,2.238281e-03,1.441406e-03,9.109375e-04...
    5.562500e-04,3.679688e-04,2.242188e-04,1.587190e-04,6.165715e-05]; %snr 30-50
BER_A_32_32_16=[5.367734e-01,5.365961e-01,5.371383e-01,5.367992e-01,5.366883e-01,5.369648e-01,5.371992e-01,...
    5.364031e-01,5.373281e-01,5.370996e-01,5.371652e-01];  %snr 30-50

%The expectation of OM is 8.230591e-01
BER_E_32_4_16=[9.231250e-02,5.450000e-02,2.618750e-02,8.750000e-03,1.423366e-03,1.444385e-04,4.312962e-06]; %snr 0-12
BER_A_32_4_16=[9.431250e-02,5.681250e-02,2.731250e-02,9.750000e-03,1.978479e-03,3.798734e-04,1.296476e-04];  %snr 0-12

%% Figure plotting
figure (1)
semilogy([0:2:8], BER_E_128_8_16, '-*', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:8], BER_A_128_8_16, ':*', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:12], BER_E_64_8_16, '-^', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:12], BER_A_64_8_16, ':^', 'MarkerSize',6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:16], BER_E_128_32_16, '-+', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:16], BER_A_128_32_16, ':+', 'MarkerSize',6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
title('16QAM MIMO systems');
hold on;
legend('128X8 MIMO MMSE EMI', '128X8 MIMO MMSE SENIA-1',...
    '64X8 MIMO MMSE EMI', '64X8 MIMO MMSE SENIA-1','128X32 MIMO MMSE EMI', ...
    '128X32 MIMO MMSE SENIA-1');
hold off;


figure (2)
semilogy([30:2:50], BER_E_32_32_16, '-*', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([30:2:50], BER_A_32_32_16, ':*', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:12], BER_E_32_4_16, '-^', 'MarkerSize', 6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
semilogy([0:2:12], BER_A_32_4_16, ':^', 'MarkerSize',6), xlabel('SNR (dB)'), ylabel('BER');
hold on;
title('16QAM MIMO systems');
hold on;
legend('32X32 MIMO MMSE EMI', '32X32 MIMO MMSE SENIA-1',...
    '32X4 MIMO MMSE EMI', '32X4 MIMO MMSE SENIA-1');
hold off;



