%% Membership functions of Interval Type-2 Fuzzy sets
% 
%                            exp(-(x1-15)^2*0.04)                               exp(-(x1-15)^2*0.04)      
%   LMF1(x1) =       -----------------------------     UMF1(x1) =       ---------------------------------
%                       0.05+exp(-(x1-15)^2*0.04)                          0.01+exp(-(x1-15)^2*0.04)
% 
%                               0.925                                       0.925  
%   LMF3(x1) = 0.025 -  ---------------------     UMF3(x1) = 0.025 -  ---------------------
%                        1+exp(-(x1-4.5)/8)                            1+exp(-(x1-3.5)/8)
%   
%                                                                  
%   LMF2(x1) = 1 - UMF1(x1) - UMF3(x1)            UMF2(x1) = 1 - LMF1(x1) - LMF3(x1)
%       
% 
%% 
function [LowerMF, UpperMF, mf_min, mf_max,mfsum_min, mfsum_max, mfsub_min, mfsub_max] = MF4_numeric(x1)
% x1 is a premise variable
x2 = x1;
LMF_1 = exp(-0.04*(x1-15).^2) ./ (0.05+exp(-0.04*(x1-15).^2));   
UMF_1 = exp(-0.04*(x1-15).^2) ./ (0.01+exp(-0.04*(x1-15).^2));   

LMF_2 = exp(-0.1*(x1).^2) ./ (0.04+exp(-0.1*(x1).^2));   
UMF_2 = exp(-0.1*(x1).^2) ./ (0.01+exp(-0.1*(x1).^2));

LMF_3 = exp(-(x2+8).^2/15) ./ (0.06+exp(-(x2+8).^2/15));   
UMF_3 = exp(-(x2+8).^2/15) ./ (0.01+exp(-(x2+8).^2/15));   

LMF_4 = exp(-((x2+16).^2)/22) ./ (0.05+exp(-((x2+16).^2)/22));   
UMF_4 = exp(-((x2+16).^2)/22) ./ (0.01+exp(-((x2+16).^2)/22));   




LowerMF = [LMF_1; LMF_2; LMF_3; LMF_4];
UpperMF = [UMF_1; UMF_2; UMF_3; UMF_4];

PlotMF(x1,LowerMF,UpperMF)

LMF1 = LMF_1.*LMF_3;
UMF1 = UMF_1.*UMF_3;

LMF2 = LMF_1.*LMF_4;
UMF2 = UMF_1.*UMF_4;

LMF3 = LMF_2.*LMF_3;
UMF3 = UMF_2.*UMF_3;

LMF4 = LMF_2.*LMF_4;
UMF4 = UMF_2.*UMF_4;

LowerMF = [LMF1; LMF2; LMF3; LMF4];
UpperMF = [UMF1; UMF2; UMF3; UMF4];

PlotMF(x1,LowerMF,UpperMF)
dim = size(UpperMF);
nRules = dim(1);
sumLowerMF = zeros(1,dim(2));
sumUpperMF = zeros(1,dim(2));

mf_min = min(LowerMF');
mf_max = max(UpperMF');

for i = 1:nRules
    sumLowerMF(1,:) = sumLowerMF(1,:) + LowerMF(i,:);
    sumUpperMF(1,:) = sumUpperMF(1,:) + UpperMF(i,:);
end
mfsum_min = min(sumLowerMF');
mfsum_max = max(sumUpperMF');

mfsub_min = min(UpperMF' - LowerMF');
mfsub_max = max(UpperMF' - LowerMF');

end





