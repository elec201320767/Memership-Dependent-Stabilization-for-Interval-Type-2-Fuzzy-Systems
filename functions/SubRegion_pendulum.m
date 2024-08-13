%% SubRegion For Type 2 Fuzzy Model

function [LowerMF, UpperMF, a, b, c, d] = SubRegion_pendulum(premise, Div, SysParam)
% SubRegion 
% input parameter 
%  : "premise" is premise variables of a fuzzy model
%    "Div" is, an integer larger than 0, the number of division of FOU (footprint of uncertainty)
%    "p" is the number of rules of a fuzzy model
f1 = SysParam(1,:);
f2 = SysParam(2,:);
mp = SysParam(3,:);
Mc = SysParam(4,:);

Domain = length(premise);
SubDomain = round(Domain/Div);
iter = 1;

if Domain ~= SubDomain * Div
    fprintf('\n Error Occured! :            == SubRegion ==\n')
    fprintf('                    The length of every subdomain should be same\n')
    return
end

while iter <= Div
    [LowerMF(:,:,iter), UpperMF(:,:,iter), a(:,:,iter), b(:,:,iter), c(:,:,iter), d(:,:,iter)] = MF4_pendulum(premise(SubDomain*(iter-1) + 1 : SubDomain*iter), f1, f2, mp, Mc);
    iter = iter + 1;
end
end





