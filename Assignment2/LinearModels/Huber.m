close all
load Problem1Data.mat
path(path,'../../exportfig');


%% Huber
m=length(y);
tau = 1.0293;  %this is equivalent to the value of the data noise std dev
I = eye(m,m);
e=ones(m,1);

A = [t, ones(length(y),1)];
A = [-A -I I I];

b = y;

H = [zeros(2*m+2,2*m+2)   zeros(2*m+2,m); 
     zeros(m,2*m+2)        eye(m,m)];

g = tau.*[zeros(2,1); e; e; zeros(m,1)]; 

%lower and upper bounds
L = [-Inf*ones(2,1);zeros(m,1);zeros(m,1);-Inf*ones(m,1)];

N = 100;
tic; %timing
for i=1:N
x=quadprog(H,g,[],[],A,-b,L); %solving
elapsedTime = toc/N;
end

A = [t, ones(length(y),1)];
x = x(1:2);
dataPlot;
hold on
plot(t,A*x,'r')
title('Huber estimation')
legend('data','true model','Huber model','location','southeast')
hold off
print('../img/Huber','-dpng')
stats(t,y,A,x,'../img/Huber_hist','../img/Huber_CI')