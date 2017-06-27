close all
load Problem1Data.mat
path(path,'../../exportfig');


%% L1 - norm

A = [t, ones(length(y),1)];
b = y;
f=[0; 0; ones(length(t),1)]; 


[m,n]=size(A);
I = eye(m,m);
M=[A I ; 
  -A I];
B=[b;-b];

N = 100;
tic; %timing
for i = 1:N
x=linprog(f,-M,-B);
end
elapsedTime = toc/N;

figure
hold on;
x=x(1:2);
dataPlot;
plot(t,A*x,'r')
legend('data','true model','L1 model','location','southeast')
print('../img/L1_norm','-dpng')

stats(t,y,A,x,'../img/L1_hist','../img/L1_CI');