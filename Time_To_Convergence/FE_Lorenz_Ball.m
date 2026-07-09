clear;clc;

epsilon = 1e-16;
sig = 10;
rho = 8/3;
beta = 28;
dt = 1e-5;
s0t = [6;5;9];
s0 = [1203;-1049;360];

MU = diag([10.9795971282064,14.6854347255662,0.00127537481573532]);



%% Defining Operators
dxdt = @(s) -sig*(s(1)-s(2));
dydt = @(s) rho*s(1)-s(2)-s(1)*s(3);
dzdt = @(s) s(1)*s(2)-beta*s(3);

dPosdt = @(s) [dxdt(s);dydt(s);dzdt(s)];

%% Calculate Analytic Bound on Time to Convergence

%% Calculate Discrete Bound on Time to Convergence



%%Running a Single Simulation
max_time = 100;

LorenzData = Forward_Euler(s0t,dPosdt,max_time,dt,[],[]);
LorenzPred = Forward_Euler(s0,dPosdt,max_time,dt,LorenzData,MU);
LorenzDiff = LorenzData-LorenzPred;

for col = 1:size(LorenzDiff,2)
    LorenzDiffNorm(1,col) = norm(LorenzDiff(:,col));
end
%%
timeline = (0:size(LorenzDiffNorm,2)-1)*dt;

semilogy(timeline,LorenzDiffNorm)

%% Not to be run. Used originally to calculate minimal MU values.

% % This script is meant to calculate the minimal nudging coefficients for
% % the Lorenz system when we know that the system will at some point be
% % within an absorbing ball (because the Lorenz '63 system is dissipative).
% clear; clc; yalmip clear;
% negdefconst = 1e-9;
% 
% % Define the decision variables
% mu = sdpvar(3,1);
% X = sdpvar(3,1);
% Xtrue = sdpvar(3,1);
% % [-2.3365;3.263;1.7192];
% optvar = mu;
% 
% 
% % Define related constants
% sig = 10;
% beta = 8/3;
% rho = 28;
% 
% % Here's all the YALMIP stuff
% 
% sym_matrix = [-(sig+mu(1)) (sig+rho-Xtrue(3))/2 Xtrue(2)/2;
%               (sig+rho-Xtrue(3))/2 -(1+mu(2)) 0;
%               Xtrue(2)/2 0 -(beta+mu(3))];
% dUdt = X'*sym_matrix*X
% 
% neg_dUdt_sos = -dUdt-negdefconst*dot(X,X); % Negative definiteness constraint
% 
% obj_func = dot(mu,mu);
% % Introduce the constraints
% % Hmm, not working so well: , X(1)^2 <= 2.3365^2, X(2)^2 <= 3.263^2, X(3) >= 0, X(3) <= 1.7192
% constr = [sos(neg_dUdt_sos), mu(:) >= 0, Xtrue(1)^2 <= 2.3365^2, Xtrue(2)^2 <= 3.263^2, Xtrue(3) >= 0, Xtrue(3) <= 1.7192];
% 
% options = sdpsettings('solver', 'mosek');
% 
% [sol,v,M,res] = solvesos(constr,obj_func,options,optvar);
% 
% sdisplay(value(mu))
