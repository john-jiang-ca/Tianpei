%% this routine plot the BER versus SNR(dB) figure for RSVD-GA
%% BER data
BER_RSVD_GA_32_4=[0.1425312500, 0.0804140625, 0.0206640625, 0.0022109375, 0.0001477192];
BER_RSVD_GA_32_16=[0.2870507812, 0.2520664062, 0.2172187500, 0.1861093750, 0.1514804688,...
    0.1044882812, 0.0398359375, 0.0083007812, 0.0015859375];

%% ploting figure
%% 4QAM
figure(1)
semilogy([4:2:12], BER_RSVD_GA_32_4, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on
semilogy([4:2:20], BER_RSVD_GA_32_16, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('BER');
hold on 
title('BER Versus SNR RSVD-GA'), legend('RSVD-GA-4QAM', 'RSVD-GA-16QAM');
hold off