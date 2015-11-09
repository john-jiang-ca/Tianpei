function [ symOut, MSE ] = CSVD_modified( y, H, SNRd, M, pav, C,  epsilon, tol, dataMod )
% this routine implement complex Support vector detector for MIMO system
% the global dual optimizaiton problem is divided into two separate dual
% optimizaiton problems
K=H*H';
% y=y+conj(y);
Nr=length(K(1,:));
Kr=real(K);
Ki=imag(K);
Yr=real(y);
Yi=imag(y);
LambdaR=zeros(Nr,1);
LambdaI=zeros(Nr,1);
%% real channel 
[LambdaR]=a2DSolver_modified(Yr, Kr, C, epsilon, tol);

%% imaginary channel
[LambdaI]=a2DSolver_modified(Yi, Kr, C, epsilon, tol);

%%Reconstruction
Lambda=complex(LambdaR, LambdaI);
% I=eye(Nr);
% G=H'*(I/K);
% G=(I/(H'*H))*H';
% symOut_mmse=G*y;
% symOut_hat=2*G*Kr*Lambda;
symOut_hat=H'*Lambda;
% u_c=H'*Lambda;   %conjugate of u
% v_c=H.'*Lambda;   %conjugate of v
% H_c=complex(real(H), -imag(H));
% symOut_hat=u_c+G*H_c*v_c;
MSE=norm(H*dataMod-2*Kr*Lambda);
% MSE2=norm(H*dataMod-H*symOut_mmse);
MSE3=norm(H*dataMod-H*symOut_hat);
symOut=Rectangular_QAM_slicer(symOut_hat,M, pav);
end