%% Setup
clear;clc;

x0t = [.1;.1;-.1;.1];
x0 = [.2;.2;-.2;.1];
sx0 = x0-x0t;
epsilon = 1e-16;
dt = 1e-5;

mumin = [0.823543613122324,0,0,0]; %

Operator = [-1 2 0.5 0.2;
             1 -2 1 1;
             0.5 0.5 -3 -1;
             0.5 0.5 1 -4];
MU = diag(mumin);
Opsize = size(Operator,1);

II = [1.328:0.00000001:1.328297];

EE = size(II,2);
posis = zeros(1,EE);

for HH = 1:EE
   [P_c,D_c] = eig(Operator-MU/II(HH));
    Discrete_Time = dt/log(norm(eye(Opsize)+D_c*dt))*log(epsilon/((norm(P_c)*norm((P_c^-1)*sx0))));
    if Discrete_Time>0
        posis(HH) = 1;
    else
        posis(HH) = 0;
    end
end

mindex = find(posis==0,1,"first");

MINMU = mumin(1)/II(mindex)

%%
% Calculating the exact value for mumin.

eigvmax = @(u) (max(real(eig(Operator-diag([u,0,0,0])))));
MINMUexact = fzero(eigvmax,MINMU);

