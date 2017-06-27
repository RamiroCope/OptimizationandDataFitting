function [ ratio ] = autocorrelationTest( r )
%AUTOCORRELATIONTEST This function computes the values of autocorrelation
%and trend treshold for any input of r residuals. This helps assess whether
%short sequences of consecutive residuals are correlated.
%
%   INPUT: r is the residuals vector
%   OUTPUT: the autocorrelation to trend treshold ratio
%
%NOTE: If the absolute value of the autocorrelation exceeds the trend
%treshold, a trend is likely to be present in the residuals. Therefore,
%values greater than 1 suggest that a trend is likely to be present in the
%residuals. 

m = length(r);
v = zeros(1, m-1);

for i=1:m-1
    v(1,i) = r(i)*r(i+1);
end

auto = abs(sum(v));

trend = (r'*r) / sqrt(m-1);

ratio = auto/trend;
end
