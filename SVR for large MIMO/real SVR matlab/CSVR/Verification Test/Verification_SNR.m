%% verification test for SVD
close all;
clear all;
tic
SNR=2:2:60;   %receive signal to noise ratio in dB
SNRd=10.^(SNR./10);   %receive signal to noise ratio decimal
M=30;   %number of training data samples
M_t=1000; %number of testing data samples
N=30; %the dimensions of the linear function
pav=1/N;      %the average power of the regression coefficiences
w=complex(normrnd(0, sqrt(pav/2), [N,1]), normrnd(0, sqrt(pav/2),[N,1]));   %the true regression coeffeciences
% Pt=norm(w)^(2);   %the total power of the true function 
d_hat=zeros(M,1); %the true output of the linear function in training part
d=zeros(M,1);   %the observations in training part
w_svr=zeros(M,1); %the regression estimation coefficiences from SVR
noiseV=1./SNRd;
realization=1e2;  %the realization time
epsilon=0;
C=10;
tol=1e-3;
MSE_training=zeros(length(SNR),1);  %training output MSE for SVR
MSE_svrco_training=zeros(length(SNR),1);  %training MSE for coefficiences vector for SVR
MSE_testing=zeros(length(SNR),1);   %testing output MSE for SVR
MSE_mmse_training=zeros(length(SNR),1); %training of the output of MSE for MMSE
MSE_mmseco_training=zeros(length(SNR),1); %training MSE for coefficiences vector for MMSE
MSE_mmse_testing=zeros(length(SNR),1); %testing MSE for MMSE

%% generate the file to record the simulation results
fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\Verification Test\Test Data\MSE_coefficientVector.txt', 'a');
fprintf(fid, '-------------\n');
fprintf(fid, 'this document record the data of SVR from testing\n');
fprintf(fid, 'this document compare the output of SVR and MMSE regression from the view of SNR\n');
fprintf(fid, 'this document also record the MSE comparison of the regression coefficience vector (coefficience vector) of SVR and MMSE\n');
fprintf(fid, 'this document consider the influence of the weight of regularization term to the estimation accuracy of coefficience vector\n');
fprintf(fid, 'all the result are based on %d realizations\n', realization);
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'The training data set size is %d with dimension %d\n', M, N);
fprintf(fid, 'The testing data set size is %d\n', M_t);
fprintf(fid, 'the hyperparameters of SVR are\n');
fprintf(fid, 'C %f \n', C);
fprintf(fid, 'epsilon %e\n', epsilon);
fprintf(fid, 'tolerence %e\n', tol);
    
for count=1:length(SNR)
    

      MSE1=0;  % MSE for SVR training
      MSE1_co=0; %MSE for SVR coefficience vector training
      MSE_training_mmse=0; %MSE for MMSE training
      MSE_training_mmse_co=0; % MSE for MMSE coefficiences vector training
      
      MSE2=0;  %MSE for SVR testing
      MSE_testing_mmse=0;% MSE for MMSE testing
 
for count1=1:realization

%% The traing part
x=complex(normrnd(0,sqrt(1/2), M, N), normrnd(0, sqrt(1/2), M ,N));    %input data matrix
n=complex(normrnd(0, sqrt(noiseV(count)/2), M,1), normrnd(0, sqrt(noiseV(count)/2), M,1));  %AWGN noise vector
d_hat=x*w;  %the true output of the linear function
d=d_hat+n; %the observations
w_svr=SVR_learning(d, x, SNRd(count),C ,epsilon, tol);   %get the regression estimation by SVD (in the form of regression coefficiences)
d_out_svr=x*w_svr;     % the  regression estimation from SVR
MSE1=MSE1+norm(d_out_svr-d_hat)^(2)/M;   % the mean MSE of the output of the training data of SVR  
MSE1_co=MSE1_co+norm(w-w_svr)^(2)/M;   %the mean MSE of the coefficience vector for the training data of SVR
% the MMSE regression
I=eye(N);
w_mmse=(I/(x'*x+SNRd(count)^(-1)*I))*x'*d; 
d_out_mmse=x*w_mmse;
MSE_training_mmse=MSE_training_mmse+norm(d_hat-d_out_mmse)^(2)/M;  %the mean MSE of the output of the training data of  MMSE 
MSE_training_mmse_co=MSE_training_mmse_co+norm(w-w_mmse)^(2)/M; %the mean MSE of the coefficience vector of MMSE training data 



%% the testing part

x_t=complex(normrnd(0, sqrt(1/2), M_t, N), normrnd(0, sqrt(1/2), M_t, N)); %get the test data set
n_t=complex(normrnd(0, sqrt(noiseV(count)/2), M_t, 1), normrnd(0, sqrt(noiseV(count)/2), M_t, 1));  % the noise of the observation of test data
d_hat_t=x_t*w;   %the true output of the test data 
d_t=d_hat_t+n_t;   %observation of the test data output
d_out_svr_t=x_t*w_svr;   %the output of the SVR regresssion estimation 
MSE2= MSE2+norm(d_out_svr_t-d_hat_t)^(2)/M_t;   % the mean MSE of the testing data samples

d_out_mmse_t=x_t*w_mmse;
MSE_testing_mmse=MSE_testing_mmse+norm(d_hat_t-d_out_mmse_t)^(2)/M_t;   %the mean MSE of the mmse regression



end
MSE_training(count)=MSE1/realization; 
MSE_svrco_training(count)=MSE1_co/realization; 
MSE_mmse_training(count)=MSE_training_mmse/realization;
MSE_mmseco_training(count)=MSE_training_mmse_co/realization;
MSE_testing(count)=MSE2/realization;
MSE_mmse_testing(count)=MSE_testing_mmse/realization;
end
fprintf(fid, '-----------------\n');
fprintf(fid, 'SVR\n');
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d, ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of training data set\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_training(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of testing data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_testing(count));
end
fprintf(fid, '\n');
fprintf(fid, '-----------------\n');



fprintf(fid, '-----------------\n');
fprintf(fid, 'MMSE\n');
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d, ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of training data set\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_mmse_training(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of testing data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_mmse_testing(count));
end
fprintf(fid, '\n');
fprintf(fid, '-----------------\n');



fprintf(fid, '-----------------\n');
fprintf(fid, 'SVR\n');
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d, ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of the coefficience vector\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_svrco_training(count));
end
fprintf(fid, '\n');
fprintf(fid, '-----------------\n');



fprintf(fid, '-----------------\n');
fprintf(fid, 'MMSE\n');
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d, ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'MSE of the coefficience vector\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_mmseco_training(count));
end
fprintf(fid, '\n');
fprintf(fid, '-----------------\n');




fprintf(fid, 'the program ends successfully\n');
fprintf(fid, '--------------\n');
fclose(fid);
   

%% figure plotting
% figure(1)
% title('Prediction accuraccy of output comparison between SVR and MMSE versus SNR');
% hold on 
% plot(SNR, MSE_testing, '-*', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('Prediction Risk');
% hold on 
% plot(SNR, MSE_mmse_testing, '--+', 'MarkerSize', 8), xlabel('SNR(dB)'), ylabel('Prediction Risk');
% hold on 
% legend('SVR', 'MMSE');
% hold off
% 
% figure(2)
% title('Training Error of output comparison between SVR and MMSE versus SNR');
% hold on 
% plot(SNR, MSE_training, '-*', 'MarkerSize', 5), xlabel('SNR(dB)'), ylabel('Training Error');
% hold on 
% plot(SNR, MSE_mmse_training, '--+', 'MarkerSize', 5), xlabel('SNR(dB)'), ylabel('Training Error');
% hold on 
% legend('SVR', 'MMSE');
% hold off
% 
% figure(3)
% title('Training Error comparison of coefficience vector between SVR and MMSE versus SNR');
% hold on 
% plot(SNR, MSE_svrco_training, '-*', 'MarkerSize', 5), xlabel('SNR(dB)'), ylabel('Training Error');
% hold on 
% plot(SNR, MSE_mmseco_training, '--+', 'MarkerSize', 5), xlabel('SNR(dB)'), ylabel('Training Error');
% hold on 
% legend('SVR', 'MMSE');
% hold off


toc