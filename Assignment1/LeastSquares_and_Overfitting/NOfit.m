function [ x_star, r_star, A ] = NOfit( t , y , n )

%NOFIT is a function that computes the least squares fit of the following
%model: M(x,t) = x1 +x2 sin(wt)+x3 cos(wt)+x4 sin(2wt)+x5 cos(2wt)+···
%given the inputs of t, y and n. Here w = 2pi/24 is the period.

%INPUTS: t is a column vector representing the x-axis (time) of the time 
%          series data in hours
%        y is a column vector with the observed measurements of [NOx]
%        n is the order of the fit model desired; that is, n is the number 
%          of features desired to be included in the approxmation model.
%            

%OUTPUTS: The function outputs:
%                  x_star: a column vector holding the n unknown parameters 
%                  r_star: a column vector holding the residuals 
%                  A     : matrix A is added to ease computations in
%                          OptimalN.m script. It is the A matrix of the 
%                          normal equations to solve the LSQ fit problem.


%We define matrix A with each column representing a basis
%function (eg. f1(t1)=x1, f2(t1)=x2*sin(wt)) and with rows representing
%the different inputs (t1, t2,...,t24) for each basis function.


w = 2*pi/24;  W = (1:length(t))'*w;
N = n-1;
A = zeros(length(t),N);
close all;
scatter(t,y);


for ii=1:2:N
    A(:,ii)=sin(W*(ceil(ii/2))); 
end

for ii=2:2:N
    A(:,ii)=cos(W*(ii/2)); 
end


A=[ones(length(t),1) A];
%We solve x_star by using the normal equation
x_star = linsolve((A'*A), (A'*y));
r_star = y - A*x_star;

hold on
plot(t,A*x_star)
end

