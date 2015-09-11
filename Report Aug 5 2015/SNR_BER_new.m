%this is for ploting the figure of SNR-BER curve
%the output of this code generate a eps file 
% sim8_4_0_18=[0.22589062154, 0.15894874930, 0.08640749753, 0.03090812452, 0.00636250013, 0.00079124997, 0.00005481876,0.00000279009,0.00000014225951464822, 0.00000000451437990952];
% sim8_16_0_18=[0.25537124276, 0.22414547205, 0.19270686805, 0.15720969439, 0.10538093746, 0.04198125005, 0.00701593747, 0.00052656251,0.00004961727, 0.00000235999];
% sim8_64_0_22=[0.30989542603, 0.28304228187, 0.25555500388, 0.22751520574, 0.19991666079, 0.16876062751, 0.12789249420, 0.07114583254, 0.02219000086, 0.00312604173, 0.00018895834,0.00000633392  ];
% sim16_4_0_14=[0.23751124740, 0.17347000539, 0.09674874693, 0.03350968659, 0.00618375000, 0.00060968747, 0.00002066603,0.00000026777  ];
% sim16_16_0_18=[0.25537124276, 0.22414547205, 0.19270686805, 0.15720969439, 0.10538093746, 0.04198125005, 0.00701593747, 0.00052656251,0.00002031656,0.00000020876338590830];
% 0.25537124276 0.22414547205 0.19270686805 0.15720969439, 0.075001, 0.019427, 0.002040, 0.000093,
% sim32_4_0_14=[0.24545359612, 0.18672686815, 0.11097452790, 0.04044124857, 0.00811937544, 0.00093781250,0.00005755865,0.00000097851];
% sim48_4_0_14=[0.25442603230, 0.19931979477, 0.12637500465, 0.05224249884, 0.01231614593, 0.00185479166, 0.00016166667,0.00000549426];
% sim_100_40_old=[0.26819999999999999396, 0.19164999999999998703, 0.11859999999999999709,...
%     0.06090000000000000274, 0.02457499999999999962, 0.00645000000000000000 ,0.00121277617675312195,0.00005600000000000000,0.00000121703063996339...
%     0.00000000607241023727 ]; %the old one with one coordinate updated
sim_100_40_new=[0.29618624999999998426, 0.21398074999999999735, 0.13460074999999999124, 0.06905224999999999558, ...
    0.02625699999999999909, 0.00652400000000000028, 0.00085649999999999995, 0.00004672285869138618,0.00000078341626839089]; %the new one with 1Dsearching
%2Dsolver
% sim_100_100_old=[0.55520000000000002682, 0.51502999999999998781, 0.45769999999999999574, 0.40034999999999998366,...
%     0.33395999999999997909, 0.26013999999999998236, 0.19559000000000001385];
sim_100_100_new=[0.55873399999999995291, 0.51256619999999997184, 0.45856750000000001677, 0.39660560000000000258,...
    0.32833479999999998222, 0.25939909999999999357, 0.19414680000000000826, 0.14058170000000000388,0.10032209999999999739]; %the new 100 times 100 system
fig1=figure;
Xsize=8.5;Ysize=11;
% Xleft=(8-Xsize)/2;Ytop=(12-Ysize)/2;
Xleft=0;
Ytop=0;
set(fig1,'Units','inches');
set(fig1,'position', [Xleft,Ytop,Xsize,Ysize]);
set(fig1,'PaperUnits','inches');
set(fig1,'PaperPosition',[Xleft,Ytop,Xsize,Ysize]);
set(gca,'fontsize',15);
figure (1)
% fig1=semilogy(0:2:18,sim8_4_0_18,'-*','MarkerSize',15);xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:18,sim8_16_0_18,'->','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% fig1=semilogy(0:2:22,sim8_64_0_22,'-d','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:14,sim16_4_0_14,'-o','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:18,sim16_16_0_18,'-s','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
% fig1=semilogy(0:2:14,sim32_4_0_14,'-p','MarkerSize',15),xlabel('SNR per bit (dB)'),ylabel('BER');
% hold on;
fig1=semilogy(6:2:22,sim_100_40_new,'-+','MarkerSize',15);xlabel('SNR (dB)'),ylabel('SER'),
hold on;
fig1=semilogy(6:2:22,sim_100_100_new,'->','MarkerSize',15);xlabel('SNR (dB)'),ylabel('SER'),
hold on;
% fig1=semilogy(6:2:24,sim_100_40_old,'-*','MarkerSize',15),xlabel('SNR (dB)'),ylabel('SER'),
% legend('8X8 4QAM','8X8 16QAM','8X8 64QAM','16X16 4QAM','16X16 16QAM','32X32 4QAM','48X48 4QAM');
legend('100X40 new', '100X100 old');
hold off;
print('-depsc','-tiff','-r300','BER_curves');
% hold off;

