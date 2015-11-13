%% this routine verify the 2-D  update rule with considering noise term in dual 
% objective function

%% parameter settings
K11=2;
K11_hat=3;
K12=3;
K22=3;
K22_hat=4;
m1=5;
m2=7;
lambda1=4;
lambda2=6;


S=K11_hat*K22_hat-K12^(2);
lambda2_new=(K11_hat*K22-K12^(2))*lambda2+lambda1*K12*(K11_hat-K11)+m2*K11_hat-m1*K12;
lambda2_new=lambda2_new/S;


lambda1_new=(K22_hat*K11-K12^(2))*lambda1+lambda2*K12*(K22_hat-K22)+m1*K22_hat-m2*K12;
lambda1_new=lambda1_new/S;

%%verification of lambda1_new an lambda2_new
lambda1_true=lambda1*K11+m1-(lambda2_new-lambda2)*K12;
lambda1_true=lambda1_true/K11_hat;
lambda2_true=lambda2*K22+m2-(lambda1_new-lambda1)*K12;
lambda2_true=lambda2_true/K22_hat;
