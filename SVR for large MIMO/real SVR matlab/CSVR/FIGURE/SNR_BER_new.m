SER_MMSE_OSIC_32_4=[0.2972062500, 0.2279687500, 0.1424250000, 0.0705062500, 0.0286875000, 0.0089000000, 0.0020500000, 0.0003006077 ];
SER_MMSE_OSIC_32_16=[0.7440250000 ,0.7111250000 ,0.6748750000 ,0.6352875000 ,0.5856500000, 0.5062062500, 0.3750312500, 0.2198000000, 0.0967562500, 0.0360125000 ];
SER_MMSE_OSIC_64_4=[0.3385843750, 0.2814031250 ,0.2082000000, 0.1259750000, 0.0598437500, 0.0233000000, 0.0060375000, 0.0010937500 ];
SER_MMSE_OSIC_64_16=[0.7628968750, 0.7342406250 ,0.7039031250 ,0.6676812500, 0.6317593750, 0.5807312500 ,0.5181187500 ,0.4029062500 ];
SER_RSVD_OMIC_32_4=[0.1529875000 ,0.0724187500, 0.0253937500, 0.0111187500, 0.0047312500 ,0.0020437500, 0.0020062500 ,0.0014250000 ];                                                                                                                      
SER_RSVD_OMIC_32_16=[0.6707875000, 0.6202625000 ,0.5760937500, 0.5362062500, 0.4985500000, 0.4672812500, 0.4484625000, 0.4391375000, 0.4282062500, 0.4167125000 ];
SER_RSVD_OMIC_64_4=[0.1571250000, 0.0588906250, 0.0090093750, 0.0008875000, 0.0001451598, 0.0000446016 ];
SER_RSVD_OMIC_64_16=[0.6781468750, 0.6220937500, 0.5678343750, 0.5115250000, 0.4639031250 ,0.4216093750 ,0.3954843750, 0.3695625000 ];
fig1=figure;
% Xsize=8.5;Ysize=11;
% % Xleft=(8-Xsize)/2;Ytop=(12-Ysize)/2;
% Xleft=0;
% Ytop=0;
% set(fig1,'Units','inches');
% set(fig1,'position', [Xleft,Ytop,Xsize,Ysize]);
% set(fig1,'PaperUnits','inches');
% set(fig1,'PaperPosition',[Xleft,Ytop,Xsize,Ysize]);
% set(gca,'fontsize',15);
figure (1)
fig1=semilogy([10:2:24],SER_MMSE_OSIC_32_4,'-*');xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:28],SER_MMSE_OSIC_32_16,'-d','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:24],SER_MMSE_OSIC_64_4,'-o','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:24],SER_MMSE_OSIC_64_16,'->','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;

fig1=semilogy([10:2:24],SER_RSVD_OMIC_32_4,':*','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:28],SER_RSVD_OMIC_32_16,':d','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:20],SER_RSVD_OMIC_64_4,':o','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
fig1=semilogy([10:2:24],SER_RSVD_OMIC_64_16,':>','MarkerSize',15);xlabel('SNR(dB)'),ylabel('SER');
hold on;
% fig1=semilogy(0:2:18,sim8_16_0_18,'->','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% fig1=semilogy(0:2:22,sim8_64_0_22,'-d','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:14,sim16_4_0_14,'-o','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:18,sim16_16_0_18,'-s','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:14,sim32_4_0_14,'-p','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(6:2:22,sim_100_40_new,'-+','MarkerSize',15);xlabel('SNR (dB)'),ylabel('SER'),
% hold on;
% fig1=semilogy(6:2:22,sim_100_100_new,'->','MarkerSize',15);xlabel('SNR (dB)'),ylabel('SER'),
% hold on;
% fig1=semilogy(6:2:24,sim_100_40_old,'-*','MarkerSize',15),xlabel('SNR (dB)'),ylabel('SER'),
% legend('8X8 4QAM','8X8 16QAM','8X8 64QAM','16X16 4QAM','16X16 16QAM','32X32 4QAM','48X48 4QAM');
legend('MMSE-OSIC-32-4', 'MMSE-OSIC-32-16', 'MMSE-OSIC-64-4', 'MMSE-OSIC-64-16',...
    'RSVD-OMIC-32-4', 'RSVD-OMIC-32-16', 'RSVD-OMIC-64-4', 'RSVD-OMIC-64-16');
hold off;
% print('-depsc','-tiff','-r300','BER_curves');
% hold off;

