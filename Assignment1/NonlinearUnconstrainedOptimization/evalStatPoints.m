X = [3, 2; -3.77931, -3.28319; -2.80512, 3.13131; 
    3.58443, -1.84813; -3.07303, -0.081353; -0.127961, -1.95371; 
    0.0866775, 2.88425; 3.38515, 0.0738519; -0.270845, -0.923039];

file = fopen('../Latex/evalStatPoints.tex','w');

for i = 1:length(X)
    [~,g,H] = himmelblau(X(i,:));
    E = eig(H);
    fprintf(file,['\\paragraph{Stationary point \\#%d:} $(x_1,x_2) = ' ...
        '(%.2f,%.2f)$\n'],i,X(i,1),X(i,2));
    fprintf(file,'\\begin{gather*}\n');
    fprintf(file,['\\nabla f(x) = \\begin{pmatrix} %.2f \\\\ %.2f ' ...
        '\\end{pmatrix}, \\quad \n'],g(1),g(2));
    fprintf(file,['\\nabla^2 f(x) = \\begin{pmatrix} %.2f & %.2f \\\\ ' ... 
        '%.2f & %.2f \\end{pmatrix} \\\\ \n'],H(1,1),H(1,2),H(2,1),H(2,2));
    fprintf(file,'\\text{Eigenvalues of Hessian: } %.2f, %.2f',E(1),E(2));
    fprintf(file,'\\end{gather*}\n');
    if sum(E>0) == 2
        fprintf(file,['The Hessian is positive definite since both ' ...
            'eigenvalues are positive. This means that the stationary ' ...
            'point $(x_1,x_2) = (%.2f,%.2f)$ is a local minimmum'], ...
            X(i,1),X(i,2));
    elseif sum(E>0) == 1
        fprintf(file,['The Hessian is indefinite since the eigenvalues '...
            'have opposite sign. This means that the stationary point '...
            '$(x_1,x_2) = (%.2f,%.2f)$ is a saddle point'],X(i,1),X(i,2));
    elseif sum(E>0) == 0
        fprintf(file,['The Hessian is negative definte since both ' ...
            'eigenvalues are negative. This means that the stationary ' ...
            'point $(x_1,x_2) = (%.2f,%.2f)$ is a local maximum'], ...
            X(i,1),X(i,2));
    end
end

fclose(file);