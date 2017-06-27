function [yhat,J] = estimateY(theta,data)
t = data(:,1);
n = 1; np = length(theta);
z0 = [10 0 0];
[T,Z] = ode45(@ModelAndSensitivity,t,z0,[],theta,n,np);
yhat = Z(:,1);
Sp = Z(:,2:3);
J = Sp;