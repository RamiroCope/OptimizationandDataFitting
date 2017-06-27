function phiContours(x,y)
% Compute the function values
x1 = linspace(0,1,1000);
x2 = linspace(0,5,1000);
[Theta1,Theta2]=meshgrid(x1,x2);

phi = zeros(size(Theta1));
for i = 1:size(Theta1,1)
    for j = 1:size(Theta1,2)
        phi(i,j) = 1/2*sum((y-Theta1(i,j)*x./(Theta2(i,j) + x)).^2);
    end
end

v = linspace(0,1^(1/3),20).^3;
contourf(Theta1,Theta2,phi,v,'linewidth',2);


xlabel('\theta_1','Fontsize',14)
ylabel('\theta_2','Fontsize',14)
set(gca,'fontsize',14);
set(gcf, 'Color', 'w');
colorbar