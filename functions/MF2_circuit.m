%% Membership functions of Interval Type-2 Fuzzy sets
% 
%                               0.925                                       0.925  
%   LMF1(x1) = 0.95 -  ---------------------     UMF1(x1) = 0.95 -  ---------------------
%                       1+exp(-(x1+4.5)/8)                            1+exp(-(x1+3.5)/8)
% 
%                               0.925                                       0.925  
%   LMF3(x1) = 0.025 -  ---------------------     UMF3(x1) = 0.025 -  ---------------------
%                        1+exp(-(x1-4.5)/8)                            1+exp(-(x1-3.5)/8)
%   
%                                                                  
%   LMF2(x1) = 1 - UMF1(x1) - UMF3(x1)            UMF2(x1) = 1 - LMF1(x1) - LMF3(x1)
%       
% 
%% Cart-Pendulum MFs with 4 rules
function [LowerMF, UpperMF, mf_min, mf_max, mfsum_min, mfsum_max, mfsub_min, mfsub_max] = MF2_circuit(x1)
% x1 is a premise variable
 
C = 15;
L = 1000;
R = 10;
b = [0.1 0.2];
v_max = max(x1);

g_min = -max(b)/C*v_max^2
g_max = 0;
 
LMF_1 = (g_max-(-min(b)/C*x1.^2))/(g_max-g_min);
UMF_1 = (g_max-(-max(b)/C*x1.^2))/(g_max-g_min);

LMF_2 = ((-max(b)/C*x1.^2)-g_min)/(g_max-g_min);
UMF_2 = ((-min(b)/C*x1.^2)-g_min)/(g_max-g_min);

plot(x1,LMF_1); hold on
plot(x1,UMF_1); hold on
plot(x1,LMF_2); hold on
plot(x1,UMF_2); hold on

LowerMF(1,:) = LMF_1;
LowerMF(2,:) = LMF_2;
UpperMF(1,:) = UMF_1;
UpperMF(2,:) = UMF_2;

mf_min = min(LowerMF');
mf_max = max(UpperMF');

dim = size(UpperMF);
nRules = dim(1);
sumLowerMF = zeros(1,dim(2));
sumUpperMF = zeros(1,dim(2));
for i = 1:nRules
    sumLowerMF(1,:) = sumLowerMF(1,:) + LowerMF(i,:);
    sumUpperMF(1,:) = sumUpperMF(1,:) + UpperMF(i,:);
end
mfsum_min = min(sumLowerMF');
mfsum_max = max(sumUpperMF');

mfsub_min = min(UpperMF' - LowerMF');
mfsub_max = max(UpperMF' - LowerMF');

end





