cd 'C:\Users\Preston Chen\Documents\GitHub\Tianpei\FCSD-CUDA project\parallel_acceleration_of_FCSD_based_on_cuda'
%this is for ploting the figure of SNR-BER curve
sim8_4_0_18=[0.22589062154, 0.15894874930, 0.08640749753, 0.03090812452, 0.00636250013, 0.00079124997, 0.00005481876,0.00000279009];
sim8_16_0_18=[0.25537124276, 0.22414547205, 0.19270686805, 0.15720969439, 0.10538093746, 0.04198125005, 0.00701593747, 0.00052656251,0.00004961727, 0.00000235999];
sim8_64_0_22=[0.30989542603, 0.28304228187, 0.25555500388, 0.22751520574, 0.19991666079, 0.16876062751, 0.12789249420, 0.07114583254, 0.02219000086, 0.00312604173, 0.00018895834,0.00000633392  ];
sim16_4_0_14=[0.23751124740, 0.17347000539, 0.09674874693, 0.03350968659, 0.00618375000, 0.00060968747, 0.00002066603,0.00000026777  ];
sim16_16_0_18=[0.25537124276, 0.22414547205, 0.19270686805, 0.15720969439, 0.10538093746, 0.04198125005, 0.00701593747, 0.00052656251,0.00002031656,0.00000020876338590830];
% 0.25537124276 0.22414547205 0.19270686805 0.15720969439, 0.075001, 0.019427, 0.002040, 0.000093,
sim32_4_0_14=[0.24545359612, 0.18672686815, 0.11097452790, 0.04044124857, 0.00811937544, 0.00093781250,0.00005755865,0.00000097851];
sim48_4_0_14=[0.25442603230, 0.19931979477, 0.12637500465, 0.05224249884, 0.01231614593, 0.00185479166, 0.00016166667,0.00000549426];

fig1=figure;
Xsize=8.5;Ysize=11;
% Xleft=(8-Xsize)/2;Ytop=(12-Ysize)/2;
Xleft=0;
Ytop=0;
set(fig1,'Units','inches');
set(fig1,'position', [Xleft,Ytop,Xsize,Ysize]);
set(fig1,'PaperUnits','inches');
set(fig1,'PaperPosition',[Xleft,Ytop,Xsize,Ysize]);
set(gca,'fontsize',20);
set(gca,'FontName', 'Arial');
figure (1)
fig1=semilogy(0:2:14,sim8_4_0_18,'-k*','MarkerSize',15);xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:18,sim8_16_0_18,'-k>','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
fig1=semilogy(0:2:22,sim8_64_0_22,'-kd','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:14,sim16_4_0_14,'-ko','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:18,sim16_16_0_18,'-ks','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:14,sim32_4_0_14,'-kp','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
hold on;
fig1=semilogy(0:2:14,sim48_4_0_14,'-k+','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER'),
legend({'8X8 4QAM','8X8 16QAM','8X8 64QAM','16X16 4QAM','16X16 16QAM','32X32 4QAM','48X48 4QAM'},'FontSize',80,'FontWeight','bold');
legend('Location','SouthWest');
% po=get(gca,'Position');
% set(gca,'Position',[po(1)-0.1,po(2),po(3),po(4)]);
hold on;
print('-depsc','-tiff','-r300','BER_curves');
hold off;

