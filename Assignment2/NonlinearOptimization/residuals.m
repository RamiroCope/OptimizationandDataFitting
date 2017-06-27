function [r,J] = residuals(theta,data)
[yhat,J] = estimateY(theta,data);
r = data(:,2) - yhat;
