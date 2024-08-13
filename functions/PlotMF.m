%% Plot membership funciotns
%% 
function PlotMF(x,LMF,UMF)
% x1 is a premise variable
temp = size(UMF);
p = temp(1);    % p is # of rules

LinS = {'-','--',':','-.'};
figure()
for i = 1 : p
    plot(x,UMF(i,:),'k','linestyle',LinS{i});hold on
    plot(x,LMF(i,:),'k','linestyle',LinS{i});hold on
end
xlim([min(x) max((x))]); ylim([0 1]);

% legend('$\theta_1(x_1(t))$ Boundary','$\theta_2(x_1(t))$ Boundary','$\theta_3(x_1(t))$ Boundary','Interpreter','Latex')
% legend('Upper and lower bound for $\theta_1$ ','Upper and lower bound for $\theta_2$','Upper and lower bound for $\theta_3$','Interpreter','Latex')
xlabel('$x_1(t)$','Interpreter','latex')
 
end
% memberhship function(=grade of membership) plot
% figure(1)
% for i = 1:1
% plot(x,UpperMF(1,:),'k');hold on
% plot(x,UpperMF(2,:),'-.','Color','k');hold on
% plot(x,UpperMF(3,:),'--','Color','k');hold on
% plot(x,LowerMF(1,:),'k');hold on
% plot(x,LowerMF(2,:),'-.','Color','k');hold on
% plot(x,LowerMF(3,:),'--','Color','k');hold on
% xlim([-10 10]); ylim([0 1]);
% legend('$\theta_1(x_1(t))$ Boundary','$\theta_2(x_1(t))$ Boundary','$\theta_3(x_1(t))$ Boundary','Interpreter','Latex')
% xlabel('$x_1(t)$','Interpreter','latex')
% end


