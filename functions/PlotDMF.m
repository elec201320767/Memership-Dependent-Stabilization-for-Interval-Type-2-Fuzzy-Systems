%% Plot membership funciotns
%% 
function [VerticesDL, VerticesDU, rho] = PlotDMF(x,step,LMF,UMF) 
% x1 is a premise variable
Nrules = size(UMF);
p = Nrules(1);    % p is # of rules
DiffUMF = diff(UMF,1,2)/step;
DiffLMF = diff(LMF,1,2)/step;
redundancy = 1;
LastDiff = length(DiffUMF);
for i = 1 : p
    DiffUMF(i,LastDiff + redundancy) = DiffUMF(i,LastDiff);
    DiffLMF(i,LastDiff + redundancy) = DiffLMF(i,LastDiff);
end 
VerticesDL = ExtremeBox([DiffLMF; DiffLMF])';
VerticesDU = ExtremeBox([DiffUMF; DiffUMF])';

rhoLlow = min(DiffLMF');
rhoLupp = max(DiffLMF');
rhoUlow = min(DiffUMF');
rhoUupp = max(DiffUMF');
rho = [rhoLlow;rhoLupp;rhoUlow;rhoUupp];

LinS = {'-','--',':','-.'};
figure()
for i = 1 : p
    plot(x,DiffUMF(i,:),'k','linestyle',LinS{i});hold on
    plot(x,DiffLMF(i,:),'k','linestyle',LinS{i});hold on
end
xlim([min(x) max((x))]);  

legend('$\dot\theta_1(x_1(t))$ Boundary','$\dot\theta_2(x_1(t))$ Boundary','$\dot\theta_3(x_1(t))$ Boundary','Interpreter','Latex')
xlabel('$x_1(t)$','Interpreter','latex')
 
end
 


