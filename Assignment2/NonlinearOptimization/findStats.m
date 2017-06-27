function [Cov, ci1, ci2] = findStats(x,y,theta,J,alpha)

m = length(x); n = length(theta);

t = tinv(alpha/2,m-n);
Var = 1/(m-n)*sum((y-theta(1)*x./(theta(2)+x)).^2);

H = J'*J;
C = diag(H^(-1));
Conf = t*sqrt(Var)*sqrt(C);
Cov = Var*H^(-1);

ci1 = sort([Conf(1) -Conf(1)] + theta(1));
ci2 = sort([Conf(2) -Conf(2)] + theta(2));
