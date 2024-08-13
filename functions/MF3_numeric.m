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
%% 
function [LowerMF, UpperMF, mf_min, mf_max,mfsum_min, mfsum_max, mfsub_min, mfsub_max] = MF3_numeric(x1)
% x1 is a premise variable
LMF_1 = 0.95 - (0.925./(1+exp(-((x1+4.5)/8))));   
UMF_1 = 0.95 - (0.925./(1+exp(-((x1+3.5)/8))));    

LMF_3 = 0.025 + (0.925./(1+exp(-((x1-4.5)/8))));
UMF_3 = 0.025 + (0.925./(1+exp(-((x1-3.5)/8))));

LMF_2 = 1 - (UMF_1 + UMF_3);
UMF_2 = 1 - (LMF_1 + LMF_3);

LowerMF = [LMF_1; LMF_2; LMF_3];
UpperMF = [UMF_1; UMF_2; UMF_3];
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





