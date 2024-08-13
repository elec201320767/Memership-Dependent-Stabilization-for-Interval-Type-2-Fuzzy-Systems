%% Sublinear Models
% 
%  A1 = [2.78 -5.63;    A2 = [0.2  -3.22;       A3 = [-avar   -6.63;
%         0.01 0.33]         0.35 0.12]                0.45 0.15] 
% 
% B1 = [2;              B2 = [8;                B3 = [-bvar+6;
%      -1]                    0]                         -1]
% 
%%
function [A, B] = SubSys3_numeric(avar,bvar)
n = 2; % # of states
m = 1; % # of input
p = 3; % # of rules

A1 = [2.78 -5.63; 0.01 0.33];
A2 = [0.2  -3.22; 0.35 0.12];
B1 = [2;    -1];
B2 = [8;     0];
A3 = [-avar   -6.63; 0.45 0.15];
B3 = [-bvar+6; -1];

 A = zeros(n,n,p); B = zeros(n,m,p);
 
 A(:,:,1) = A1;A(:,:,2) = A2;A(:,:,3) = A3;
 B(:,:,1) = B1;B(:,:,2) = B2;B(:,:,3) = B3;
   
end





