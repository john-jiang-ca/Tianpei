%% verification test for SVD
close all;
clear all;
tic
SNR=2:2:40;   %receive signal to noise ratio in dB
SNRd=10.^(SNR./10);   %receive signal to noise ratio decimal
M=30;   %number of training data samples
M_t=1000; %number of testing data samples
N=30; %the number of features
pav=1/N;      %the average power of the regression coefficiences
w=complex(normrnd(0, sqrt(pav/2), [N,1]), normrnd(0, sqrt(pav/2),[N,1]));   %the true regression coeffeciences
% Pt=norm(w)^(2);   %the total power of the true function 
d_hat=zeros(M,1); %the true output of the linear function in training part
d=zeros(M,1);   %the observations in training part
w_svr=zeros(M,1); %the regression estimation coefficiences from SVR
noiseV=1./SNRd;
realization=1e2;  %the realization time
epsilon=1e-7;
C=1;
tol=1e-3;
MSE_training=zeros(length(SNR),1);  %training MSE for SVR
MSE_testing=zeros(length(SNR),1);   %testing MSE for SVR
MSE_mmse_training=zeros(length(SNR),1); %training MSE for MMSE
MSE_mmse_testing=zeros(length(SNR),1); %testing MSE for MMSE
fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\Verification Test\Test Data\MSE_sparse_data_setting.txt', 'a');
fprintf(fid, '-------------\n');
fprintf(fid, 'this document record the data of SVR from testing\n');
fprintf(fid, 'this document compare the SVR and MMSE from the view of data sparseness');
fprintf(fid, 'all the result are based on %d realizations', realization);
fprintf(fid, 'the SNR(dB) are\n');
for count=1:length(SNR)
    fprintf(fid, '%d ', SNR(count));
end
fprintf(fid, '\n');
fprintf(fid, 'The training data set size is %d with %d feature\n', M, N);
fprintf(fid, 'The testing data set size is %d\n', M_t);
fprintf(fid, 'the hyperparameters of SVR are\n');
fprintf(fid, 'C %f \n', C);
fprintf(fid, 'epsilon %0.10f\n', epsilon);
fprintf(fid, 'tolerence %0.10f\n', tol);

for count=1:length(SNR)
      MSE1=0;
      MSE_training_mmse=0;
      MSE2=0;
      MSE_testing_mmse=0;
for count1=1:realization

%% The traing part
x=complex(normrnd(0,sqrt(1/2), M, N), normrnd(0, sqrt(1/2), M ,N));    %input data matrix
n=complex(normrnd(0, sqrt(noiseV(count)/2), M,1), normrnd(0, sqrt(noiseV(count)/2), M,1));  %AWGN noise vector
d_hat=x*w;  %the true output of the linear function
d=d_hat+n; %the observations
w_svr=SVR_learning(d, x, SNRd(count),C ,epsilon, tol);   %get the regression estimation by SVD (in the form of regression coefficiences)
d_out=x*w_svr;     % the  regression estimation from SVR
MSE1=MSE1+norm(d_out-d_hat)^(2)/M;   % the mean MSE of the training data   


% the MMSE regression
I=eye(N);
w_mmse=(I/(x'*x+SNRd(count)^(-1)*I))*x'*d; 
d_out_mmse=x*w_mmse;
MSE_training_mmse=MSE_training_mmse+norm(d_hat-d_out_mmse)^(2)/M;



%% the testing part

x_t=complex(normrnd(0, sqrt(pav/2), M_t, N), normrnd(0, sqrt(pav/2), M_t, N)); %get the test data set
n_t=complex(normrnd(0, sqrt(noiseV(count)/2), M_t, 1), normrnd(0, sqrt(noiseV(count)/2), M_t, 1));  % the noise of the observation of test data
d_hat_t=x_t*w;   %the true output of the test data 
d_t=d_hat_t+n_t;   %observation of the test data output
d_out_t=x_t*w_svr;   %the output of the SVR regresssion estimation 
MSE2= MSE2+norm(d_out_t-d_hat_t)^(2)/M_t;   % the mean MSE of the testing data samples

d_out_mmse_t=x_t*w_mmse;
MSE_testing_mmse=MSE_testing_mmse+norm(d_hat_t-d_out_mmse_t)^(2)/M_t;   %the mean MSE of the mmse regression



end
MSE_training(count)=MSE1/realization; 
MSE_mmse_training(count)=MSE_training_mmse/realization;
MSE_testing(count)=MSE2/realization;
MSE_mmse_testing(count)=MSE_testing_mmse/realization;
end
fprintf(fid, '-----------------\n');
fprintf(fid, 'SVR\n');
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





fprintf(fid, 'the program ends successfully\n');
fprintf(fid, '--------------\n');
fclose(fid);


%% figure plotting
% figure(1)





toc

