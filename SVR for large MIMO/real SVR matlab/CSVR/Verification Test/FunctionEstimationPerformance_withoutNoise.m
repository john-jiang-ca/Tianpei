%% Test of SVR and MMSE under the condition that without additive noise
close all;
clear all;
tic
% SNR=2:2:50;   %receive signal to noise ratio in dB
% SNRd=10.^(SNR./10);   %receive signal to noise ratio decimal
M=30;   %number of training data samples
M_t=100; %number of testing data samples
N=30; %the dimensions of the linear function
pav=1/N;      %the average power of the regression coefficiences
% Pt=norm(w)^(2);   %the total power of the true function 
d_hat=zeros(M,1); %the true output of the linear function in training part
d=zeros(M,1);   %the observations in training part
w_svr=zeros(M,1); %the regression estimation coefficiences from SVR
% noiseV=1./SNRd;
realization=1e3;  %the realization time
epsilon=0;
C=1e15;
tol=1e-3;
SNR=0;
MSE_training_svrV=zeros(length(SNR),1);  %training output MSE for SVR
MSE_training_svr_coV=zeros(length(SNR),1);  %training MSE for coefficiences vector for SVR
MSE_testing_svrV=zeros(length(SNR),1);   %testing output MSE for SVR
MSE_training_mmseV=zeros(length(SNR),1); %training of the output of MSE for MMSE
MSE_training_mmse_coV=zeros(length(SNR),1); %training MSE for coefficiences vector for MMSE
MSE_testing_mmseV=zeros(length(SNR),1); %testing MSE for MMSE

%% generate the file to record the simulation results
fid=fopen('F:\GitHub\Tianpei\SVR for large MIMO\real SVR matlab\CSVR\Verification Test\Test Data\FunctionEstimation_withoutNoise.txt', 'a');
fprintf(fid, '=======================================\n');
fprintf(fid, 'this document records the comparison of the accuracy of function estimation between SVR and MMSE without additive noise\n');
% fprintf(fid, 'this document compare the output of SVR and MMSE regression from the view of SNR\n');
% fprintf(fid, 'this document also record the MSE comparison of the regression coefficience vector (coefficience vector) of SVR and MMSE\n');
% fprintf(fid, 'this document consider the influence of the weight of regularization term to the estimation accuracy of coefficience vector\n');
fprintf(fid, 'SYSTEM CONFIGURATION\n');
fprintf(fid, '*********************\n');
fprintf(fid, 'all the result are based on %d realizations\n', realization);
% fprintf(fid, 'the SNR(dB) are\n');
% for count=1:length(SNR)
%     fprintf(fid, '%d ', SNR(count));
% end
% fprintf(fid, '\n');
fprintf(fid, 'The training data set size is %d with dimension %d\n', M, N);
fprintf(fid, 'The testing data set size is %d\n', M_t);
fprintf(fid, 'the hyperparameters of SVR are\n');
fprintf(fid, 'C %f \n', C);
fprintf(fid, 'epsilon %e\n', epsilon);
fprintf(fid, 'tolerence %e\n', tol);
fprintf(fid, '*********************\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
% for count=1:length(SNR)
    

      MSE_training_svr=0;  % MSE for SVR training output
      MSE_training_svr_co=0; %MSE for SVR coefficience vector training
      MSE_training_mmse=0; %MSE for MMSE training output
      MSE_training_mmse_co=0; % MSE for MMSE coefficiences vector training
      
      MSE_testing_svr=0;  %MSE for SVR testing output
      MSE_testing_mmse=0;% MSE for MMSE testing output
 
for count1=1:realization

%% The traing part
w=complex(normrnd(0, sqrt(pav/2), [N,1]), normrnd(0, sqrt(pav/2),[N,1]));   %the true regression coeffeciences   
x=complex(normrnd(0,sqrt(1/2), M, N), normrnd(0, sqrt(1/2), M ,N));    %input data matrix
% n=complex(normrnd(0, sqrt(noiseV(count)/2), M,1), normrnd(0, sqrt(noiseV(count)/2), M,1));  %AWGN noise vector
d_hat=x*w;  %the true output of the linear function
d=d_hat; %the observations
w_svr=SVR_learning(d, x, SNR,C ,epsilon, tol);   %get the regression estimation by SVD (in the form of regression coefficiences)
d_out_svr=x*w_svr;     % the  regression estimation from SVR
MSE_training_svr=MSE_training_svr+norm(d_out_svr-d_hat)^(2)/M;   % the mean MSE of the output of the training data of SVR  
MSE_training_svr_co=MSE_training_svr_co+norm(w-w_svr)^(2)/M;   %the mean MSE of the coefficience vector for the training data of SVR
% the MMSE regression
I=eye(N);
% w_mmse=(I/(x'*x+SNRd(count)^(-1)*I))*x'*d;
w_mmse=I/(x'*x)*x'*d;
d_out_mmse=x*w_mmse;
MSE_training_mmse=MSE_training_mmse+norm(d_hat-d_out_mmse)^(2)/M;  %the mean MSE of the output of the training data of  MMSE 
MSE_training_mmse_co=MSE_training_mmse_co+norm(w-w_mmse)^(2)/M; %the mean MSE of the coefficience vector of MMSE training data 



%% the testing part

x_t=complex(normrnd(0, sqrt(1/2), M_t, N), normrnd(0, sqrt(1/2), M_t, N)); %get the test data set
% n_t=complex(normrnd(0, sqrt(noiseV(count)/2), M_t, 1), normrnd(0, sqrt(noiseV(count)/2), M_t, 1));  % the noise of the observation of test data
d_hat_t=x_t*w;   %the true output of the test data 
% d_t=d_hat_t+n_t;   %observation of the test data output
d_out_svr_t=x_t*w_svr;   %the output of the SVR regresssion estimation 
MSE_testing_svr= MSE_testing_svr+norm(d_out_svr_t-d_hat_t)^(2)/M_t;   % MSE of the prediction output SVR

d_out_mmse_t=x_t*w_mmse;
MSE_testing_mmse=MSE_testing_mmse+norm(d_hat_t-d_out_mmse_t)^(2)/M_t;   %MSE of the prediction output of MMSE



end
MSE_training_svrV=MSE_training_svr/realization; 
MSE_training_svr_coV=MSE_training_svr_co/realization; 
MSE_training_mmseV=MSE_training_mmse/realization;
MSE_training_mmse_coV=MSE_training_mmse_co/realization;
MSE_testing_svrV=MSE_testing_svr/realization;
MSE_testing_mmseV=MSE_testing_mmse/realization;
% end
fprintf(fid, 'SIMULATION OUTPUT\n');
fprintf(fid, '*********************\n');
fprintf(fid, 'Results for SVR\n');
fprintf(fid, '---------------\n');
fprintf(fid, 'output MSE of training data set\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_training_svrV(count));
end
fprintf(fid, '\n');
fprintf(fid, 'output MSE of testing data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_testing_svrV(count));
end
fprintf(fid, '\n');
fprintf(fid, 'coefficient vector MSE of training data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_training_svr_coV(count));
end
fprintf(fid, '\n');
fprintf(fid, '---------------\n');
fprintf(fid, '\n');
fprintf(fid, '\n');


fprintf(fid, 'Results for MMSE\n');
fprintf(fid, '---------------\n');
fprintf(fid, 'output MSE of training data set\n');
for count=1:length(SNR)
fprintf(fid, '%0.10f, ', MSE_training_mmseV(count));
end
fprintf(fid, '\n');
fprintf(fid, 'output MSE of testing data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_testing_mmseV(count));
end
fprintf(fid, '\n');
fprintf(fid, 'coefficient vector MSE of training data set\n');
for count=1:length(SNR)
    fprintf(fid, '%0.10f, ', MSE_training_mmse_coV(count));
end
fprintf(fid, '\n');
fprintf(fid, '---------------\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '*********************\n');
fprintf(fid,'\n');
fprintf(fid,'\n');






fprintf(fid, 'the program ends successfully\n');
fprintf(fid, '=======================================\n');
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