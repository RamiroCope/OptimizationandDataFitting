function [ z ] = runTest( r )
%RUNTEST runTest is a function used to assess the randomness of signs (+
%and -) in residuals. 
%   INPUT: A column vector of residuals wished to be analyzed
%   OUTPUT: A Z-score, that if less than 1.96, the residuals' signs can be
%           considered random at a 5% significance level. 

m=length(r);
pvector = zeros(1,m); nvector = zeros(1,m);
runs = 1;

for i=1:m
    pvector(i) = r(i) >= 0;
    nvector(i) = r(i) <= 0;
end

npos = sum(pvector); nneg = sum(nvector);

for ii=1:m-1
    if ~(pvector(ii) == pvector(ii+1))
        runs = runs +1;
    end
end

mu_runs = (2 * npos * nneg / m) + 1;
std_runs = sqrt(((mu_runs - 1) * (mu_runs - 2) / (m - 1)));

z = abs((runs - mu_runs)) / std_runs;
end


