clear; close all; clc;
path(path,'../../exportfig');
path(path,'..');

x = [2.0000 2.0000 0.6670 0.6670 0.4000 0.4000 0.2860 0.2860 0.2220 0.2200 0.2000 0.2000];
y = [0.0615 0.0527 0.0334 0.0334 0.0138 0.0258 0.0129 0.0183 0.0083 0.0169 0.0129 0.0087];
x = x'; y = y';
alpha = 0.05;
reps = 100;
timeL = zeros(1,reps); timeN = timeL; timeNJ = timeN;
theta0 = [0.0615 0.6670];

%% 2.1 Linear estimation of theta
% Plot raw data
figure('Position', [100, 0, 1000,1000]);
plot(x,y,'.','markersize',30);
xlabel('Concentration of substrate (x)','Fontsize',14)
ylabel('Reaction rate (y)','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
export_fig '../img/2-1raw.png'

% Estimate model, after transforming to Linear problem
a = 1./x; b = 1./y;
A = [a, ones(length(b),1)];
H = 2*A'*A;
g = -2*A'*b;
for i = 1:reps
    T = tic;
    [est,~,~,outputL] = quadprog(H,g,[],[],[],[],[0 0],[inf inf]); %solving
    timeL(i) = toc(T);
end

% Plot reciprocal of data with model
figure('Position', [100, 0, 1000,1000]);
plot(a,b,'.','markersize',30);
hold on
plot(a,A*est,'b');
xlabel('1/x','Fontsize',14)
ylabel('1/y','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
export_fig '../img/2-1reciprocal.png'

% Transform parameters back
thetaL(1) = 1/est(2);
thetaL(2) = est(1)*thetaL(1);
phiL = 1/2*sum((y-thetaL(1)*x./(thetaL(2) + x)).^2);

%% 2.2 Non linear estimation

%Without using jacobian
for i = 1:reps
    T = tic;
    [thetaN, resnormN, ~, ~, outputN, ~, jacobianN] = lsqnonlin(@michaelisMenten,theta0,[0 0],[inf inf],[],[x y]);
    timeN(i) = toc(T);
end

%Including Jacobian
opts = optimoptions('lsqnonlin','Jacobian','on');
for i = 1:reps
    T = tic;
    [thetaNJ, resnormNJ, ~, ~, outputNJ, ~, jacobianNJ] = lsqnonlin(@michaelisMenten,theta0,[0 0],[inf inf],opts,[x y]);
    timeNJ(i) = toc(T);
end

fig = figure('Position', [100, 0, 1000,1000]);
phiContours(x,y)
export_fig '../img/2-2phi.png'
hold on
h1 = plot(thetaL(1),thetaL(2),'b.','markersize',40);
h2 = plot(thetaN(1),thetaN(2),'r.','markersize',40);
h3 = plot(thetaNJ(1),thetaNJ(2),'g.','markersize',40);
legend([h1 h2 h3],'Linear','Nonlinear','Nonlinear with Jacobian','location','southeast')
export_fig '../img/2-2phiPoints.png'

% Plot results
figure('Position', [100, 0, 1000,1000]);
plot(x,y,'.','markersize',30);
hold on
xp = linspace(min(x),max(x));
h1 = plot(xp,thetaL(1)*xp./(thetaL(2)+xp),'b','LineWidth',3);
h2 = plot(xp,thetaN(1)*xp./(thetaN(2)+xp),'r','LineWidth',3);
h3 = plot(xp,thetaNJ(1)*xp./(thetaNJ(2)+xp),'g','LineWidth',3);
xlabel('x','Fontsize',14)
ylabel('y','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
legend([h1 h2 h3],'Linear','Nonlinear','Nonlinear with Jacobian','location','southeast')
export_fig '../img/2-2est.png'

phiN = resnormN/2; phiNJ = resnormNJ/2;
fprintf('Linear solution: \n \t phi: %.2e \n \t Average time over %d repetitions: %.4f s \n',phiL,reps,mean(timeL));
fprintf('Without Jacobian: \n \t phi: %.2e \n \t Average time over %d repetitions: %.4f s \n \t iter: %d \n \t func calls: %d \n',phiN,reps,mean(timeN),outputN.iterations,outputN.funcCount);
fprintf('With Jacobian: \n \t phi: %.2e \n \t Average time over %d repetitions: %.4f s \n \t iter: %d \n \t func calls: %d \n',phiNJ,reps,mean(timeNJ),outputNJ.iterations,outputNJ.funcCount);

% Find stats and save tables
jacobianL = [-x./(thetaL(2)+x) thetaL(1)*x./(thetaL(2)+x).^2];

[covL, ci1L, ci2L] = findStats(x,y,thetaL,jacobianL,alpha);
[covN, ci1N, ci2N] = findStats(x,y,thetaN,full(jacobianN),alpha);
[covNJ, ci1NJ, ci2NJ] = findStats(x,y,thetaNJ,full(jacobianNJ),alpha);
saveStatsTable('2linear','$\theta_1$','$\theta_2$',thetaL(1),thetaL(2),ci1L,ci2L,covL);
saveStatsTable('2nonlinear','$\theta_1$','$\theta_2$',thetaN(1),thetaN(2),ci1N,ci2N,covN);
saveStatsTable('2nonlinearJacobian','$\theta_1$','$\theta_2$',thetaNJ(1),thetaNJ(2),ci1NJ,ci2NJ,covNJ);

% Save table comparing the methods
file = fopen('../tables/2comparison.tex','w');
fprintf(file,'\\begin{tabular}{l|cccc} \\hline \\hline \n');
fprintf(file,'Algorithm & $\\phi$ & Time [ms] & Iterations & Funciton calls  \\\\ \\hline');
fprintf(file,'\\texttt{quadprog} & %.3e & %.1f & %d & - \\\\ \n',phiL,mean(timeL)*1000,outputL.iterations);
fprintf(file,'\\texttt{nonlinsq} & %.3e & %.1f & %d & %d \\\\ \n',phiN,mean(timeN)*1000,outputN.iterations,outputN.funcCount);
fprintf(file,'\\texttt{nonlinsq} with Jacobian & %.3e & %.1f & %d & %d \\\\ \n',phiNJ,mean(timeNJ)*1000,outputNJ.iterations,outputNJ.funcCount);
fprintf(file,'\\hline \\hline \n');
fprintf(file,'\\end{tabular} \n');
fclose(file);

%% 2.3
clear; close all;
load MMBatchData.mat
opts = optimoptions('lsqnonlin','Jacobian','on');
figure('Position', [100, 0, 1000,1000]);
plot(data(:,1),data(:,2),'.','markersize',30);
hold on
theta01 = [0.4 0.4]; theta02 = [0.1 0.8]; theta03 = [0.09 0.7];

[thetaODE1, ~, ~, ~, ~, ~, jacobian1] = lsqnonlin(@residuals,theta01,[],[],opts,data);
yhat1 = estimateY(thetaODE1(1,:),data);
h1 = plot(data(:,1),yhat1,'LineWidth',3);
[thetaODE2, ~, ~, ~, ~, ~, jacobian2] = lsqnonlin(@residuals,theta02,[],[],opts,data);
yhat2 = estimateY(thetaODE2(1,:),data);
h2 = plot(data(:,1),yhat2,'LineWidth',3);
[thetaODE3, ~, ~, ~, ~, ~, jacobian3] = lsqnonlin(@residuals,theta03,[],[],opts,data);
yhat3 = estimateY(thetaODE3(1,:),data);
h3 = plot(data(:,1),yhat3,'LineWidth',3);

legend([h1 h2 h3],'\theta_0 = [0.4 0.4]','\theta_0 = [0.1 0.8]','\theta_0 = [0.09 0.7]','location','northeast')

xlabel('t','Fontsize',14)
ylabel('x','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
export_fig '../img/2-3.png'

fprintf('trust-region-reflective \n')
fprintf('[0.4 0.4]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE1(1),thetaODE1(2),1/2*sum((data(:,2) - yhat1).^2));
fprintf('[0.1 0.8]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE2(1),thetaODE2(2),1/2*sum((data(:,2) - yhat2).^2));
fprintf('[0.09 0.7]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE3(1),thetaODE3(2),1/2*sum((data(:,2) - yhat3).^2));

fig = figure('Position', [100, 0, 1000,1000]);
phiContours2(data)
hold on
h1 = plot(thetaODE1(1),thetaODE1(2),'b.','markersize',40);
h2 = plot(thetaODE2(1),thetaODE2(2),'r.','markersize',40);
h3 = plot(thetaODE3(1),thetaODE3(2),'g.','markersize',40);
x = linspace(0.07,0.17,100); a = 68; b = -5.3;

plot(x,a*x + b,'y','LineWidth',2)
legend([h1 h2 h3],'\theta_0 = [0.4 0.4]','\theta_0 = [0.1 0.8]','\theta_0 = [0.09 0.7]','location','southeast')
axis([0.07 0.17 0.05 5])
export_fig '../img/2-3Phi.png'


% Get stats and save tables
alpha = 0.05;
[cov1, ci11, ci21] = findStats2(data,thetaODE1,full(jacobian1),alpha);
[cov2, ci12, ci22] = findStats2(data,thetaODE2,full(jacobian2),alpha);
[cov3, ci13, ci23] = findStats2(data,thetaODE3,full(jacobian3),alpha);
saveStatsTable('2-3(0.4 0.4)','$\theta_1$','$\theta_2$',thetaODE1(1),thetaODE1(2),ci11,ci21,cov1);
saveStatsTable('2-3(0.1 0.8)','$\theta_1$','$\theta_2$',thetaODE2(1),thetaODE2(2),ci12,ci22,cov2);
saveStatsTable('2-3(0.09 0.7)','$\theta_1$','$\theta_2$',thetaODE3(1),thetaODE3(2),ci13,ci23,cov3);


% Try with different methods
opts = optimoptions('lsqnonlin','Jacobian','on','Algorithm','levenberg-marquardt');

thetaODE1 = lsqnonlin(@residuals,theta01,[],[],opts,data);
yhat1 = estimateY(thetaODE1(1,:),data);
thetaODE2 = lsqnonlin(@residuals,theta02,[],[],opts,data);
yhat2 = estimateY(thetaODE2(1,:),data);
thetaODE3 = lsqnonlin(@residuals,theta03,[],[],opts,data);
yhat3 = estimateY(thetaODE3(1,:),data);

fprintf('LM \n')
fprintf('[0.4 0.4]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE1(1),thetaODE1(2),1/2*sum((data(:,2) - yhat1).^2));
fprintf('[0.1 0.8]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE2(1),thetaODE2(2),1/2*sum((data(:,2) - yhat2).^2));
fprintf('[0.09 0.7]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE3(1),thetaODE3(2),1/2*sum((data(:,2) - yhat3).^2));

opts = optimoptions('lsqnonlin','Jacobian','on','StepTolerance',1e-15,'MaxFunctionEvaluations',10000,'MaxIterations',10000);

thetaODE1 = lsqnonlin(@residuals,theta01,[],[],opts,data);
yhat1 = estimateY(thetaODE1(1,:),data);
thetaODE2 = lsqnonlin(@residuals,theta02,[],[],opts,data);
yhat2 = estimateY(thetaODE2(1,:),data);
thetaODE3 = lsqnonlin(@residuals,theta03,[],[],opts,data);
yhat3 = estimateY(thetaODE3(1,:),data);

fprintf('trust-region-reflective, StepTolerance = 1e-15, MaxFunctionEvaluations = 100000, MaxIterations = 10000 \n')
fprintf('[0.4 0.4]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE1(1),thetaODE1(2),1/2*sum((data(:,2) - yhat1).^2));
fprintf('[0.1 0.8]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE2(1),thetaODE2(2),1/2*sum((data(:,2) - yhat2).^2));
fprintf('[0.09 0.7]: \n \t theta1: %.2e \n \t theta2: %.2e \n \t phi %.2e \n',thetaODE3(1),thetaODE3(2),1/2*sum((data(:,2) - yhat3).^2));
