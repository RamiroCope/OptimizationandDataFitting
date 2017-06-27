%Backtracking determines the maximum step length that gives a sufficient
%decrease for a specific search direction.

%Input parameters:
%fun=Handle to the function that needs to be minimized
%x=Current x
%f=Current f(x)
%p=Step vector (search direction)
%a=Initial step length

%Ouput parameters:
%xn: new x. xn=x+a*p
%fn: new f. f(xn)
%gn: new g. g(xn)
%evals: number of function evaluations used

function [xn,fn,gn,evals] = backtracking(fun,x,f,p,a)
p = p/norm(p);
phi = @(a)(fun(x'+a*p)); %computes the value of the function for x+a*p
                         %where a=steplenght and p=direction
stop=false;
evals=0;
iter=0;
maxIter=100;
c=0.9; %Expected decrease coefficient.
t=0.7; %step length reduction coefficient.

% The step length is being reduced in each iteration while the
% Goldstein-Armijo condition it is not satisfied.
while ~stop && iter<maxIter 
    evals = evals + 1;
    if phi(a)> c*f %Goldstein-Armijo condition. Determines the expected decrease.
    a=t*a;
    else stop=true;
    end
    iter=iter+1;
end
xn=x+a*p'; %computes the new x
[fn, gn] = fun(xn);
evals = evals + 1;