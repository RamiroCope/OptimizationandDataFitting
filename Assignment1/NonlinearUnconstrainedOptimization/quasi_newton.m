clear; clc; close all;
path(path,'../../immoptibox');
path(path,'../../exportfig');

x0 = [1,1; 1,-1; -1,0]; %List of starting points
N = size(x0,1); %Number of starting points
tol = 1e-12; %tolerance of solution
maxIter = 10000; %maximum number of iterations before terminating
iter1 = zeros(1,N); iter2 = zeros(1,N); %Vectors for number of iterations
evals1 = zeros(1,N); evals2 = zeros(1,N); %Vectors for number of func evals
rep = 100; %Number of repetitions
T1 = zeros(rep,N); T2 = zeros(rep,N); %Vectors for convergence time

e1 = cell(1,N); e2 = cell(1,N);

a=5; %Backtracking step length

%Contour plot
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours

for i = 1:N
    plot(x0(i,1),x0(i,2),'k.','markersize',40)
    
    for j = 1:rep
        iter1(i) = 0; evals1(i) = 0;
        iter2(i) = 0; evals2(i) = 0;
        
        %Backtracking
        T = tic;
        X1 = x0(i,:);
        [F1,g1] = himmelblau(X1);
        B = eye(size(x0,2));
        
        while norm(g1,'inf')>tol && iter1(i) < maxIter
            p=-B\g1;
            [xn, fn, gn, evals] = backtracking(@himmelblau,X1(end,:),F1(end),p,a);
            s = xn - X1(end,:);
            y = gn - g1;
            Bs = B*s';
            B = B - (Bs/dot(s,Bs))*Bs' + (y/dot(y,s))*y';
            X1= [X1; xn]; F1 = [F1; fn];
            g1 = gn;
            iter1(i) = iter1(i) + 1;
            evals1(i) = evals1(i) + evals;
        end
        T1(j,i) = toc(T);
        
        %line search optibox
        T = tic;
        X2 = x0(i,:);
        [F2,g2] = himmelblau(X2);
        B = eye(size(x0,2));
        
        while norm(g2,'inf')>tol && iter2(i) < maxIter
            p=-B\g2;
            [xn, fn, gn, info] = linesearch(@himmelblau,X2(end,:),F2(end),g2,p);
            s = xn' - X2(end,:);
            y = gn - g2;
            Bs = B*s';
            B = B - (Bs/dot(s,Bs))*Bs' + (y/dot(y,s))*y';
            X2 = [X2; xn']; F2 = [F2; fn];
            g2 = gn;
            iter2(i) = iter2(i) + 1;
            evals2(i) = evals2(i) + info(3);
        end
        T2(j,i) = toc(T);
    end
    h1 = plot(X1(:,1),X1(:,2),'k-o','linewidth',2);
    h2 = plot(X2(:,1),X2(:,2),'r-o','linewidth',2);
    
    e1{i} = sqrt(sum((X1 - X1(end,:)).^2,2));
    e2{i} = sqrt(sum((X2 - X2(end,:)).^2,2));
end
T1 = mean(T1); T2 = mean(T2);

legend([h1,h2],'backtracking','soft line search','location','southeast')

export_fig '../img/quasiNewton.png'

errorplot('../img/quasiNewtonError.png',e1,e2,iter1,iter2,N,maxIter);

saveTable('quasiNewton.tex',x0,iter1,iter2,evals1,evals2,T1,T2)
