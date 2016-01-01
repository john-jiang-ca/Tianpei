%% Test for High performance matrix inversion approximation by exploting Neumann series expansion and 
% Newton iteration
% this is the test for the approximation performance of the refined seriers
% expansion matrix inversion scheme
% the approximation performance is evaluated by approximation residual 
% based on making average over maxTime times realizations
Nr=128;
Nt=32;
maxTime=1e3;   %the maximum realizations
Eom=1;    %the expectation of orthogonality measure
for count=2:Nt 
    Eom=Eom*(Nr-count+1)/Nr;
end
% ARV=zeros(Nt);   %the residual between experimental and theoretical approximation 
ARE=zeros(Nt); %the average approximation residual
ER=zeros(Nt); %exact average inversion residual 
%redisual
for count=1:maxTime
H=complex(normrnd(0, sqrt(1/2), Nr, Nt), normrnd(0, sqrt(1/2), Nr, Nt));
W=H'*H;    %Wishart matrix
I=eye(Nt);   %identity matrix
W_E=I/(W);   %exact inversion of W
ER=ER+I-W_E*W;   %exact matrix inverse residual 
[W_A,W_inter1]=SENIA(W);  %approximation of the matrix inversion of W
ARE=ARE+I-W_A*W; %exact approximation residual
ART=W_inter1^(4);%theoretical approximation residual
% ARV=ARV+ARE-ART;    %the accumulated error of the approximation residual 
end
% ARV=ARV./maxTime;
ARE=ARE./maxTime;
ER=ER./maxTime;


