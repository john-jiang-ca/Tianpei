X=[0.46,0.09,0.07;1.56,3.2,3.57;31,74,78;38,112,122;611,2836,3033];    %this is the operation time by second
Y=[8,16,20,25,32];   %this is the array size
figure(1)
plot(1:5,X(:,1),'-^');
hold on
figure(1)
plot(1:5,X(:,2),'-o');
hold on
figure(1)
plot(1:5,X(:,3),'-*'); legend('GPU','CPU-modigaliani','CPU-monet');
hold off
bar(X);


