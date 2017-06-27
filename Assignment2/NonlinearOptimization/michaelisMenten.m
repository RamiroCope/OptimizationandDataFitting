function [r, J] = michaelisMenten(theta,data)
x = data(:,1); y = data(:,2);
yhat = theta(1)*x./(theta(2) + x);
r = y - yhat;

J = zeros(length(x),length(theta));

if nargout > 1
   J(:,1) = -x./(theta(2)+x);
   J(:,2) = theta(1)*x./(theta(2)+x).^2;
end