function muplot(fileName,sup1,sup2,iter1,iter2,N,maxIter)

figure('Position', [100, 0, 1000,1000]);
for i = 1:N
    subplot(3,1,i);
    if iter1(i)<maxIter; plot(sup1{i},'k.-','markersize',20); end
    hold on
    if iter2(i)<maxIter; plot(sup2{i},'r.-','markersize',20); end
    xlabel('Iteration','Fontsize',14)
    ylabel('\mu','Fontsize',14)
    legend('backtracking','soft line search')
    set(gca,'fontsize',14);
end
set(gcf, 'Color', 'w');
export_fig(fileName)