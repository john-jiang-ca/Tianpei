%% RSVD-OMIC-local search
BER_RSVD_OMIC_32_4=[0.1700156250, 0.1309687500, 0.0880468750, 0.0442031250, 0.0176875000, 0.0043281250, 0.0011818910];
BER_RSVD_OMIC_64_4=[0.1752578125, 0.1358515625, 0.0948906250, 0.0489609375, 0.0162187500, 0.0029765625, 0.0002882744];


figure(1)
semilogy([4:2:16], BER_RSVD_OMIC_32_4, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on
semilogy([4:2:16], BER_RSVD_OMIC_64_4, '-o', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('SER');
hold on 
title('BER Versus SNR 4QAM'), legend('RSVD-OMIC-32', 'RSVD-OMIC-64');
hold off