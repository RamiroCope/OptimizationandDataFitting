function phiContours2(data)
% Compute the function values
x1 = linspace(0.07,0.17,100);
x2 = linspace(0.05,5,100);
[Theta1,Theta2]=meshgrid(x1,x2);
y = data(:,2);

phi = zeros(size(Theta1));
for i = 1:size(Theta1,1)
    for j = 1:size(Theta1,2)
        yhat = estimateY([Theta1(i,j) Theta2(i,j)],data);
        phi(i,j) = 1/2*sum((y-yhat).^2);
    end
end

v = linspace(3.9^(1/3),5.5^(1/3),25).^3;
contourf(Theta1,Theta2,phi,v,'linewidth',1);

xlabel('\theta_1','Fontsize',14)
ylabel('\theta_2','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
colorbar