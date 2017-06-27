%We first create the x-axis vector (time in hours) and the y-axis vector
%(observations):
t = 1:24; t=t';
y = [91.22 28.04 22.91 26.65 42.96 101.05 202.36 328.02 364.12 299.23 ...
    238.00 227.49 218.03 223.62 238.75 271.26 267.72 251.32 230.04 ... 
    206.69 170.77 131.67 143.85 157.57]';

%The A matrices and X* vectors corresponding to each n can be found by
%getting subsets from the A matrix and X* vector for n=23. That is, the
%X* vector corresponding to n=3, for instance, is the first 3 elements of
%the X* vector of n=23. 

%Therefore, we can simply calculate the A matrix and the X* vector for
%n=23, and extract subsets from them to find the matrices and vectors for 
%lower values of n:

[X23 R23 A23] = NOfit(t, y, 23);


%We cannot do the same with the residual vectors (R) since their values
%differ with each different order of n.

%Therefore, we create matrix R with column vectors storing the residuals 
%for each different odd value of n from n=3 to 23.
R = zeros(24,11);      
for i=1:11              
    n= 2*i + 1;
    R(:,i)=NOfitR(t,y,n);
end
R1 = NOfitR(t,y,1);
R = [R1 R];  %We add residuals of n=1

%We multiply x_star(s) with matrices A to obtain the approximations to the
%observations (y_hat) for each value of n. We store them in Y-hat matrix:
Y_hat = zeros(24,11);          
for i=1:11
    j = 2*i +1;
    Y_hat(:, i) = A23(:,1:j) * X23(1:j,:);
end

%%
% Part1_2: TESTING THE SOFTWARE
%
%List the found solution x* for n=3 and plot the least squares fit M(x*,t)
%together with the data points.


scatter(t,y);
hold on; 
title('LSQ Model Fit');
xlabel('Time (hrs)');
ylabel('NO Concentration (microg/m3)');
plot(t,Y_hat(:,1), 'k', 'LineWidth', 2);     % N = 3
plot(t,Y_hat(:,2), 'k');                     % N = 5
plot(t,Y_hat(:,3), 'b', 'LineWidth', 2);     % N = 7
plot(t,Y_hat(:,4), 'r', 'LineWidth', 2);     % N = 9
plot(t,Y_hat(:,5), 'c', 'LineWidth', 2);     % N = 11
plot(t,Y_hat(:,6), 'k');                     % N = 13
plot(t,Y_hat(:,7), 'k');                     % N = 15
plot(t,Y_hat(:,8), 'k');                     % N = 17
plot(t,Y_hat(:,9), 'k');                     % N = 19
plot(t,Y_hat(:,10), 'k');                    % N = 21
plot(t,Y_hat(:,11), 'k');                    % N = 23
hold off; 

%X* solution for N = 3: [186.8058;-44.9364;-93.4298]
 X_star_N3 = X23(1:3);

%%
% Part1_3: OPTIMAL ORDER OF FIT
%
%We calculate the z-scores for randomness of signs and the autocorrelation/
%trend treshold ratio:
zscores = zeros(12,1);
for i = 1:12
    zscores(i,1) = runTest(R(:,i));
end

autocorr = zeros(12,1);
for i = 1:12
    autocorr(i,1) = autocorrelationTest(R(:,i));
end

%%Plotting residuals for different values of n with their corresponding
%%z-score (for randomness of signs) and autocorrelation/trend treshold
%%ratio.
subplot(6,2,1);                                   %N=1
scatter(t,R(:,1));
title('N = 1     Z = 4.156     A/t = 4.137');


subplot(6,2,2);                                   %N=3
scatter(t,R(:,2));
title('N = 3     Z = 3.329     A/t = 3.525');

subplot(6,2,3);
scatter(t,R(:,3));                                %N=5
title('N = 5     Z = 2.488     A/t = 2.507');

subplot(6,2,4);
scatter(t,R(:,4));                                %N=7
title('N = 7     Z = 0.417     A/t = 0.792');

subplot(6,2,5);
scatter(t,R(:,5));                                %N=9
title('N = 9     Z = 0.385     A/t = 0.499');

subplot(6,2,6);
scatter(t,R(:,6));                                %N=11
title('N = 11     Z = 0.035     A/t = 0.348');

subplot(6,2,7);
scatter(t,R(:,7));                                %N=13
title('N = 13     Z = 2.568     A/t = 2.504');

subplot(6,2,8);
scatter(t,R(:,8));                                %N=15
title('N = 15     Z = 2.293     A/t = 2.510');

subplot(6,2,9);
scatter(t,R(:,9));                                %N=17
title('N = 17     Z = 2.137     A/t = 3.804');

subplot(6,2,10);
scatter(t,R(:,10));                               %N=19
title('N = 19     Z = 3.757     A/t = 4.608');

subplot(6,2,11);
scatter(t,R(:,11));                               %N=21
title('N = 21     Z = 3.819     A/t = 4.716');

subplot(6,2,12);
scatter(t,R(:,12));                               %N=23
title('N = 23     Z = 4.598     A/t = 4.596');

%%
%Part1_3b: SCALED-RESIDUAL NORM 
%
%We calculate the 2 norm of the residuals and the scaled residuals for 
%different values of n:
Rnorm = zeros(12,1);
for i=1:12
    Rnorm(i,1) = norm(R(:,i));    %Residual Norm as a function of n
end

N = 1:2:23; N = N';
M = 24; MN = M - N;           %MN is the denominator (m - n) of S* equation

SRnorm = zeros(12,1);
for i=1:12
    SRnorm(i,1) = Rnorm(i,1)       /  sqrt(MN(i,1));  %Scaled Residual Norm
end

%We plot Residuals and Scaled Residuals as function of N
SResidualsPlot = figure;
subplot(2,1,1);
scatter(N, Rnorm);
title('Residual norm as function of N');

subplot(2,1,2);
scatter(N, SRnorm);
title('Scaled-Residual norm as function of N');

%%
% Part1_4: ESTIMATING STD DEV OF SOLUTION COEFFICIENTS (x*)
%
%We compute the std deviations for the elements of x* (for x* of n = 11):
cov_x11 = (SRnorm(6,1))^2 * (A23(:,1:11)'*A23(:,1:11))^-1;
sigma_x11 = sqrt(diag(cov_x11));





