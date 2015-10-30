%% this subroutine implement the figure plotting of RSVD-GA algorithm 

%% RSVD-GA
BER_RSVD_GA_32_4=[0.1453593750, 0.0865625000, 0.0354375000,0.0072656250,0.0014061186];


%% figure plotting
figure(1)
semilogy([4:2:12], BER_RSVD_GA_32_4, '-+', 'MarkerSize', 8);xlabel('SNR(dB)'),ylabel('BER')
hold on 
legend('CSVD-GA-32')
hold on 
title('4QAM BER Versus SNR(dB)')
hold off