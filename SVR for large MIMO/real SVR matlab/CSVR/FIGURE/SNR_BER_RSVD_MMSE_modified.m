%% routine to plot BER vs SNR figure 
% RSVD modified version:
%modification includes: correct the calculation of noiseterm (use absolute value)
%reduce clipping operation
%get feasible hyper parameter settings (C=4, epsilon=1e-7, tol=1e-3)
%apply and veritfy update rule to Theta and G
%modify stopping criteria (use absolute value and use rough adaptive stopping cirteria to reduce step size when climb to the peak)
%the data come from the file '/home/tchen44/code/Tianpei/SVR for large MIMO/real SVR matlab/CSVR/test data/BER_16QAM'
%% MMSE 
BER_MMSE_32_4=[0.1736,    0.1506,    0.1255 ,   0.1012 ,   0.0783 ,   0.0589 ,   0.0424  ,  0.0309 ,   0.0205 ,   0.0147 ,   0.0096];
BER_MMSE_32_16=[ 0.1925,    0.1856,    0.1768 ,   0.1664 ,   0.1545 ,   0.1406  ,  0.1236 ,   0.1080 ,   0.0884 ,   0.0710,    0.0530];
BER_MMSE_64_4=[ 0.2034 ,   0.1824,    0.1594  ,  0.1344 ,   0.1104 ,   0.0861 ,...
    0.0661,    0.0494,    0.0359 ,   0.0243,    0.0171 ];
BER_MMSE_64_16=[0.2006,    0.1945,    0.1881,    0.1795 ,   0.1707,    0.1592 ,   0.1462 ,   0.1324 ,   0.1156 ,   0.0976,    0.0790];
BER_MMSE_128_4=[0.2258710938, 0.2010390625, 0.1711562500, 0.1412148438, 0.1154960938,...
    0.0892734375, 0.0718203125, 0.0494765625, 0.0357968750];
BER_MMSE_128_16=[0.3429824219, 0.3219785156, 0.2968339844, 0.2692324219, 0.2457695313, 0.2177773437,...
    0.1974101562, 0.1665976563, 0.1422792969, 0.1182109375, 0.0915292969];
%% RSVD
BER_RSVD_32_4=[0.1273,    0.0825 ,   0.0505 ,   0.0216,    0.0095 ,   0.0033 ,   0.0016,    0.0006 ,   0.0004, 0.0003445312];
BER_RSVD_32_16=[ 0.1768,    0.1631,    0.1386 ,   0.1122 ,   0.0868 ,   0.0616 ,   0.0443 ,   0.0320 ,...
    0.0234, 0.0176187500, 0.0174750000, 0.0167015625];
BER_RSVD_64_4=[0.1607,    0.1168,    0.0717,    0.0361,    0.0151 ,   0.0048 ,   0.0013  ,  0.0003 ,   0.0001  , 5.5550e-05,...
    0.0000327843, 0.0000293565, 0.0000243851, 0.0000200975];
BER_RSVD_64_16=[0.1885,    0.1744 ,   0.1556 ,   0.1300 ,   0.1041 ,   0.0763 ,   0.0499 ,   0.0317 ,   0.0203,...
    0.0153398438, 0.0115371094, 0.0099433594, 0.0090000000, 0.0085937500];
BER_RSVD_128_4=[0.2174023438, 0.1646230469, 0.1107285156, 0.0635839844, 0.0295097656, 0.0100957031,...
    0.0026894531, 0.0004199219, 0.0000613226, 0.0000121096];
BER_RSVD_128_16=[0.3355605469, 0.2924082031, 0.2438554687, 0.1927636719, 0.1486835938, 0.1084316406,...
    0.0710468750, 0.0431542969, 0.0232714844, 0.0125097656, 0.0072070312];


%% ploting figure
%% 4QAM
figure(1)
semilogy([10:2:30]-3, BER_MMSE_32_4, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
semilogy([10:2:30]-3, BER_MMSE_64_4, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:26], BER_MMSE_128_4, '-+', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:28], BER_RSVD_32_4, '--*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:36], BER_RSVD_64_4, '--o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
semilogy([10:2:28], BER_RSVD_128_4, '--+', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
title('BER Versus SNR 4QAM'), legend('MMSE-32', 'MMSE-64', 'MMSE-128', 'RSVD-32', 'RSVD-64', 'RSVD-128');
hold off

%% 16QAM
figure(2)
semilogy([10:2:30]-3, BER_MMSE_32_16, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
semilogy([10:2:30]-3, BER_MMSE_64_16, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:30], BER_MMSE_128_16, '-+', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:32], BER_RSVD_32_16, '--*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
semilogy([10:2:36], BER_RSVD_64_16, '--o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
semilogy([10:2:30], BER_RSVD_128_16, '--+', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
title('BER Versus SNR 16QAM'), legend('MMSE-32', 'MMSE-64', 'MMSE-128', 'RSVD-32', 'RSVD-64', 'RSVD-128');
hold off



