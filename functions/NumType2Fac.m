%% Calculation Num of permutations


function Num = NumType2Fac(m,N)
% m : # of rules
% N : order

Num = factorial(2*m + N - 1)/(factorial(N)*factorial(2*m - 1));
end






 