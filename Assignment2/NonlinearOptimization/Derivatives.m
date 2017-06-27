function [dfdx,dfdp] = Derivatives(t,x,p)
dfdx = p(1)*p(2)/(p(2)+x).^2;
dfdp(1) = x/(p(2)+x);
dfdp(2) = p(1)*x/(p(2)+x).^2;