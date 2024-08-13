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
function [LowerMF, UpperMF, mf_min, mf_max, mfsum_min, mfsum_max, mfsub_min, mfsub_max] = MF4_pendulum(x1,f1,f2,mp,Mc)
% x1 is a premise variable
param = [1/(min(mp)+min(Mc)) 1/(max(mp)+min(Mc));
         1/(min(mp)+max(Mc)) 1/(max(mp)+max(Mc))];
L = 0.5;            % the length of the pendulum
g = 9.8;            % gravity
f1_1 = (g-param(1,2)*max(mp)*L*  0^2*    cos(x1))./(4*L/3-param(1,2)*max(mp)*L*cos(x1).^2).*(sin(x1)./x1);
f1_2 = (g-param(1,2)*max(mp)*L*  5^2*    cos(x1))./(4*L/3-param(1,2)*max(mp)*L*cos(x1).^2).*(sin(x1)./x1);
f2_1 = -param(2,2)*cos(x1)./(4*L/3-param(2,2)*max(mp)*L*cos(x1).^2);
f2_2 = -param(1,1)*cos(x1)./(4*L/3-param(1,1)*min(mp)*L*cos(x1).^2);


%%% Define Grade of membership function  
muL1_1 =  (-f1_1+max(f1))/(max(f1)-min(f1));
muL2_1 =  muL1_1;
muU3_1 =  (+f1_1-min(f1))/(max(f1)-min(f1));
muU4_1 =  muU3_1;
muU1_1 = (-f1_2+max(f1))/(max(f1)-min(f1));
muU2_1 =  muU1_1;
muL3_1 =  (+f1_2-min(f1))/(max(f1)-min(f1));
muL4_1 =  muL3_1;

muL1_2 =  (-f2_1+max(f2))/(max(f2)-min(f2));
muL3_2 =  muL1_2;
muU2_2 =  (+f2_1-min(f2))/(max(f2)-min(f2));
muU4_2 =  muU2_2;
muU1_2 =  (-f2_2+max(f2))/(max(f2)-min(f2));
muU3_2 =  muU1_2;
muL2_2 =  (+f2_2-min(f2))/(max(f2)-min(f2));
muL4_2 =  muL2_2;
%%% Define membership function 
LMF_1 = muL1_1.*muL1_2;
UMF_1 = muU1_1.*muU1_2;
LMF_2 = muL2_1.*muL2_2;
UMF_2 = muU2_1.*muU2_2;
LMF_3 = muL3_1.*muL3_2;
UMF_3 = muU3_1.*muU3_2;
LMF_4 = muL4_1.*muL4_2;
UMF_4 = muU4_1.*muU4_2;
     
LowerMF = [LMF_1; LMF_2; LMF_3; LMF_4];
UpperMF = [UMF_1; UMF_2; UMF_3; UMF_4];

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





