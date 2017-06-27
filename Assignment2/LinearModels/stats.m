function [noise_var,noise_std,x_std,x_cov] = stats(t,y,A,x,filename1,filename2)

%% HISTOGRAMS

e=y-A*x; %errors
figure
subplot(2,1,1)
histogram(e,40)

A(e>5,:) = [];
y(e>5) = [];
t(e>5) = [];
e(e>5) = []; 


subplot(2,1,2)
histogram(e,40)

set(gcf, 'Color', 'w');

print(filename1,'-dpng')
%% STATISTICAL ANALYSIS

m = length(y);    %Number of observations
n = length(x);    %Number of parameters

invAA = (A'*A)^-1;

noise_var   = (e'*e)/ (m - n);
noise_std = sqrt(noise_var); 

x_cov = noise_var * invAA; 
x_std = sqrt(diag(x_cov));

%% Confidence Intervals - PARAMETERS
alpha = 0.05; % alpha
ts = tinv(1-alpha/2, m-n); % t-student)
    
Upperx(1) = x(1) + (ts * noise_std * (invAA(1,1))^0.5);
Lowerx(1) = x(1) - (ts * noise_std * (invAA(1,1))^0.5);

Upperx(2) = x(2) + (ts * noise_std * (invAA(2,2))^0.5);
Lowerx(2) = x(2) - (ts * noise_std * (invAA(2,2))^0.5);

x1_CL = [Lowerx(1), Upperx(1)]; %For correct format for table function
x2_CL = [Lowerx(2), Upperx(2)];

%% Confidence Intervals - PREDICTIONS
yhat = A*x;

UpperPred = zeros(m,1);
LowerPred = zeros(m,1);
for i=1:m
    UpperPred(i) = yhat(i) - (ts * noise_std * ( A(i,:) * invAA * A(i,:)')^0.5);
    LowerPred(i) = yhat(i) + (ts * noise_std * ( A(i,:) * invAA * A(i,:)')^0.5);
end
% Prediction interval
% 
% UpperPred1 = zeros(m,1);
% LowerPred1 = zeros(m,1);
% for i=1:m
%     UpperPred1(i) = yhat(i) - (ts * noise_std * ( A(i,:) * invAA * A(i,:)'+1)^0.5);
%     LowerPred1(i) = yhat(i) + (ts * noise_std * ( A(i,:) * invAA * A(i,:)'+1)^0.5);
% end


%% Plot confidence intervals

figure
hold on

h1 = plot(t,A*x,'r'); % Model

h2 = plot(t,A*Lowerx','Color',[0.8,0.5,0],'LineStyle','--'); % Confidence inteval of the paramemeters
plot(t,A*Upperx','Color',[0.8,0.5,0],'LineStyle','--');

h3 = plot(t,LowerPred,'Color',[0,0.5,0],'LineStyle','--'); % Confidence interval of the predictions
plot(t,UpperPred,'Color',[0,0.5,0],'LineStyle','--');

% h4 = plot(t,LowerPred1, 'm--'); % Prediction interval
% plot(t,UpperPred1, 'm--')

set(gcf, 'Color', 'w');

e=y-A*x;
y = y(e < 10);
t = t(e < 10); 
scatter(t,y,6,'b')
legend([h1 h2 h3],'model','95% CI parameters','95% CI predictions','location','southeast')
print(filename2,'-dpng')