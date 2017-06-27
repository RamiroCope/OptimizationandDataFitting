function [f,g,H]= himmelblau(X)
%The flags are used in order to avoid unnecessary computations
x1=X(1);
x2=X(2);
f = (x1^2+x2-11)^2 + (x1+x2^2-7)^2;
g = [4*(x1^2+x2-11)*x1+2*(x1+x2.^2-7); 2*(x1^2+x2-11)+4*(x1+x2^2-7)*x2];
H = [12*x1^2+4*x2-42,       4*(x1+x2);
      4*(x1+x2),      4*x1+12*x2^2-26]; 
