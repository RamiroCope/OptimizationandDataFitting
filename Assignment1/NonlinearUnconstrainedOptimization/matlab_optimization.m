clear; clc;
x0 = [1,1; 1,-1; -1,0]; %List of starting points
N = size(x0,1); %Number of starting points
rep = 100; %Number of repetitions
opts = optimoptions('fminunc','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-12);
file = fopen('fminuncOutput.txt','w');

for i = 1:N
    T1 = zeros(rep,1);
    for j = 1:rep
        T = tic;
        [x1,fval1,exitflag1,output1]=fminunc(@himmelblau,x0(i,:),opts);
        T1(j) = toc(T);
    end
    T1 = mean(T1)*1000;
    fprintf(file,'Start: (%d,%d) Converged to: (%.2f,%.2f) Iterations: %d Evals: %d Time: %.1f Using: %s \n',x0(i,1),x0(i,2),x1(1),x1(2),output1.iterations,output1.funcCount,T1,output1.algorithm);
end
T1 = mean(T1);
fclose(file);