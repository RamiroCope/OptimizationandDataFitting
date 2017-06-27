close all
load Problem1Data.mat
path(path,'../../exportfig');

%% L-infinity - norm

A = [t, ones(length(y),1)];
b = y;
f=[0; 0; 1];


[m,n]=size(A);
e = ones(m,1);
M=[A e ; 
  -A e];
B=[b;-b];

N = 100;
tic; %timing
for i = 1:100
x=linprog(f,-M,-B);
end
elapsedTime = toc/N;

figure
hold on;
x=x(1:2);
dataPlot;
plot(t,A*x,'r')
title('L-infinity estimation')
legend('data','true model','L-infinity model','location','southeast')
print('../img/LInf_norm','-dpng')

stats(t,y,A,x,'../img/LInf_hist','../img/LInf_CI');
