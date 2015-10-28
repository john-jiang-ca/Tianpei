%% routine to plot SER vs SNR figure 
% RSVD modified version:
%modification includes: correct the calculation of noiseterm (use absolute value)
%reduce clipping operation
%get feasible hyper parameter settings (C=4, epsilon=1e-7, tol=1e-3)
%apply and veritfy update rule to Theta and G
%modify stopping criteria (use absolute value and use rough adaptive stopping cirteria to reduce step size when climb to the peak)
%the data come from the file '/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/test data/SER_16QAM'
%% MMSE 
SER_MMSE_32_4=[0.3471625000, 0.3011687500, 0.2510812500, 0.2024062500, 0.1565687500, 0.1178125000, 0.0847562500,...
    0.0618375000, 0.0410187500, 0.0293312500 ,0.0191062500];
SER_MMSE_32_16=[0.7700937500, 0.7424125000, 0.7070312500 ,0.6658000000, 0.6181000000, 0.5622250000,...
    0.4945250000, 0.4318125000, 0.3537750000 ,0.2838937500, 0.2119875000];
SER_MMSE_64_4=[0.4067812800, 0.3648468780, 0.3187862800, 0.2688406280 ,0.2208128000, 0.1722062800 ,0.1322280000, 0.0987093780,...
    0.0718893780, 0.0486878000, 0.0341812800 ];
SER_MMSE_64_16=[0.8025093750, 0.7781718750, 0.7522625000 ,0.7180750000, 0.6829937500, 0.6368843750, 0.5848312500 ,0.5295250000, 0.4625500000,...
    0.3902406250, 0.3158062500];

%% RSVD
SER_RSVD_32_4=[0.2546250000, 0.1650625000, 0.1009375000, 0.0431875000, 0.0189062500, 0.0066250000,...
    0.0031562500, 0.0011417611, 0.0007253946];
SER_RSVD_32_16=[0.7070312500, 0.6522500000, 0.5545937500, 0.4486875000, 0.3472812500,...
    0.2463750000, 0.1773125000, 0.1281875000, 0.0935625000];
SER_RSVD_64_4=[0.3213750000, 0.2335625000, 0.1433906250, 0.0721718750,0.0301406250, 0.0095312500, 0.0026093750, 0.0006688784,...
    0.0002301517, 0.0001110993];
SER_RSVD_64_16=[0.7541250000, 0.6975781250, 0.6224687500, 0.5200937500, 0.4164218750, 0.3052812500, 0.1994375000,...
    0.1268281250, 0.0812968750];



%% ploting figure
%% 4QAM
figure(1)
semilogy([10:2:30]-3, SER_MMSE_32_4, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on
semilogy([10:2:30]-3, SER_MMSE_64_4, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on 
semilogy([10:2:26], SER_RSVD_32_4, ':*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on 
semilogy([10:2:28], SER_RSVD_64_4, ':o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on
title('SER Versus SNR 4QAM'), legend('MMSE-32', 'MMSE-64', 'CSVD-32', 'CSVD-64');
hold off

%% 16QAM
%% 4QAM
figure(2)
semilogy([10:2:30]-3, SER_MMSE_32_16, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on
semilogy([10:2:30]-3, SER_MMSE_64_16, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on 
semilogy([10:2:26], SER_RSVD_32_16, ':*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on 
semilogy([10:2:26], SER_RSVD_64_16, ':o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on
title('SER Versus SNR 16QAM'), legend('MMSE-32', 'MMSE-64', 'CSVD-32', 'CSVD-64');
hold off



