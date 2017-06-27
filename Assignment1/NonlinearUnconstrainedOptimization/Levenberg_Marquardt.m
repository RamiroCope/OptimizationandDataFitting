clear; clc; close all;
path(path,'../../exportfig');

x0 = [1,1; 1,-1; -1,0]; %List of starting points
N = size(x0,1); %Number of starting points
M = size(x0,2);
tol = 1e-12; %tolerance of solution
maxIter = 10000; %maximum number of iterations before terminating
iter1 = zeros(1,N); %Vectors for number of iterations
evals1 = zeros(1,N); %Vector for number of func evals
rep = 100; %Number of repetitions
T1 = zeros(rep,N); %Vectors for convergence time

e1 = cell(1,N);

a=5; %Backtracking step length

%Contour plot
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours

for i = 1:N
    plot(x0(i,1),x0(i,2),'k.','markersize',40)
    
    for j = 1:rep
        iter1(i) = 0; evals1(i) = 0;
        T = tic;
        X=x0(i,:);
        mu=0.2;
        [F,g,H]=himmelblau(X);
        H=H+mu*eye(M);
        
        while norm(g,'inf')>tol && iter1(i)<maxIter
            [L,a]=chol(H,'lower');
            while a~=0 %adjust mu until H positive definite
                mu=2*mu;
                H=H+mu*eye(M);
                [L,a]=chol(H,'lower');
            end
            p=L'\(L\(-g)); %calculate direction H*h=-g
            xn=X(end,:)+p';
            [fn]=himmelblau(xn);
            evals1(i) = evals1(i) + 1;
            q = F(end,:)+p'*g+(1/2).*p'*H*p;
            ro=(F(end,:)-fn)/(F(end,:)-q);
            if F(end,:)-fn>0
                X=[X;xn];
                F=[F;fn];
                iter1(i)=iter1(i)+1;
                mu=mu*max(1/3,1-(2*ro-1)^3);
            else
                mu=mu*2;
            end
            [f,g,H]=himmelblau(X(end,:));
            evals1(i) = evals1(i) + 1;
            H=H+mu*eye(M);
        end
        T1(j,i) = toc(T);
    end
    h1 = plot(X(:,1),X(:,2),'k-o','linewidth',2);
    
    e1{i} = sqrt(sum((X - X(end,:)).^2,2));
end
T1 = mean(T1);
export_fig '../img/LM.png'

% Error plot
figure('Position', [100, 0, 1000,1000]);
for i = 1:N
    subplot(3,1,i);
    if iter1(i)<maxIter; semilogy(e1{i},'k.-','markersize',20); end
    xlabel('Iteration','Fontsize',14)
    ylabel('2-norm of error','Fontsize',14)
    set(gca,'fontsize',14);
end
set(gcf, 'Color', 'w');
export_fig('../img/LMError.png')

% Save table
T1 = 1000*T1;
file = fopen('../tables/LM.tex','w');
fprintf(file,'\\begin{tabular}{ccc|ccc|ccc} \\hline \\hline \n');
fprintf(file,'\\multicolumn{3}{c}{$x_0 = (%d,%d)$} & \\multicolumn{3}{c}{$x_0 = (%d,%d)$} & \\multicolumn{3}{c}{$x_0 = (%d,%d)$} \\\\ \n',x0(1,1),x0(1,2),x0(2,1),x0(2,2),x0(3,1),x0(3,2));
fprintf(file,'Iter & Evals & Time & Iter & Evals & Time & Iter & Evals & Time \\\\ \\hline \n');
fprintf(file,'%d & %d & %.2g & %d & %d & %.2g & %d & %d & %.2g \\\\ \n',iter1(1),evals1(1),T1(1),iter1(2),evals1(2),T1(2),iter1(3),evals1(3),T1(3));
fprintf(file,'\\hline \\hline \n');
fprintf(file,'\\end{tabular} \n');
fclose(file);