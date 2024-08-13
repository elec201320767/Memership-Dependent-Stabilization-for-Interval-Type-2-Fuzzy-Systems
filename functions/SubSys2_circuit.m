%% Sublinear Models


%%
function [A, B, Bw, C] = SubSys2_circuit 
n = 2; % # of states
m = 1; % # of input
p = 2; % # of rules


C = 15;
L = 1000;
R = 10;
b = [0.1 0.2];
vcMax = 3.3;
a = 0.002;

A1 = [-a/C-max(b)/C*vcMax^2 -1/C;-1/L -R/L];
A2 = [-a/C 1/C;-1/L -R/L];
B1 = [0;    1/L];
B2 = [0;    1/L];

A = zeros(n,n,p); B = zeros(n,m,p); Bw = zeros(n,m,p);
 
 A(:,:,1) = A1;A(:,:,2) = A2; 
 B(:,:,1) = B1;B(:,:,2) = B2; 
 Bw(:,:,1) = B1;Bw(:,:,2) = B2; 
 C = [1 0];
end





