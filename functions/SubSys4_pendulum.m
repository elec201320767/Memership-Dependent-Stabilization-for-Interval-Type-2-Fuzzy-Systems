%% Sublinear Models

%% Cart-Pendulum Subsystems with 4 rules
function [A, B, Bw, C] = SubSys4_pendulum(SysParam)
n = 2; % # of states
p = 1; % # of input
m = 4; % # of rules
z = 1; % # of output
f1 = SysParam(1,:);
f2 = SysParam(2,:);
mp = SysParam(3,:);
Mc = SysParam(4,:);

A1 = [0 1;min(f1) 0];   A2 = A1;
A3 = [0 1;max(f1) 0];   A4 = A3;
B1 = [0; min(f2)];      B3 = B1;
B2 = [0; max(f2)];      B4 = B2;
Bw1 = [0; 0.1]; Bw2 = [0; 0.1]; Bw3 = [0; 0.1]; Bw4 = [0; 0.1];


A = zeros(n,n,m); B = zeros(n,p,m); Bw = zeros(n,p,m); C = [z, n];
A(:,:,1) = A1;A(:,:,2) = A2;A(:,:,3) = A3; A(:,:,4) = A4;
B(:,:,1) = B1;B(:,:,2) = B2;B(:,:,3) = B3; B(:,:,4) = B4;
Bw(:,:,1) = Bw1;Bw(:,:,2) = Bw2;Bw(:,:,3) = Bw3; Bw(:,:,4) = Bw4;   
C = [1 1];
end





