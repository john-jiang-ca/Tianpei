time=1e4;
o_v=zeros(1,time);
o_v_lr=zeros(1,time);
Var_v=zeros(1,time);
Var_o=zeros(1,time);
suc=0;
 for count=1:time
[ k,W_ratio,W_LR_ratio,orthogonal,orthogonal_LR,Var,Var_o ] = channelHardening_LR(24,8 );
o_v(count)=orthogonal;
o_v_lr(count)=orthogonal_LR;
Var_v(count)=Var;
Var_o(count)=Var_o;
suc=suc+k;
 end
suc=suc/time;
o_v=round(o_v*1e3)/1e3;
o_v_lr=round(o_v_lr*1e3)/1e3;
table1=tabulate(o_v(:));
table2=tabulate(o_v_lr(:));
figure (1), plot(table1(:,1),table1(:,3)),hold on;
plot(table2(:,1),table2(:,3),'-r'),hold off;
table3=tabulate(Var_v);
table4=tabulate(Var_o);
figure (2), plot(table3(:,1),table3(:,3)),hold on;plot(table4(:,1),table4(:,3),'-r'),hold off;