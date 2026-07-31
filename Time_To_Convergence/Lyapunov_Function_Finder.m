% This script is meant to find Lyapunov Functions for particular systems of
% differential equations.

% Currently does not work.
clear; clc; yalmip clear;

%% Constants
sig = 10; beta = 8/3; rho = 28;
epsln = 1e-8;
mu = sdpvar(3,1); %[330;0;0];
MU = diag(mu);

Oper = @(s) [-sig*(s(1)-s(2)); rho*s(1)-s(2)-s(1)*s(3); s(1)*s(2)-beta*s(3)];
numvars = 3;

%% Construct Lyapunov Function Form
Xtru = sdpvar(numvars,1); % The reference system's variables (may be bounded in some way)
X = sdpvar(numvars,1); % A model of the system
dX = sdpvar(numvars,1); % Difference variables


ddXdt = Oper(Xtru+dX)-Oper(Xtru)-MU*dX;

% Monoms = monolist(dX,2,2); % Monomial of difference variables
% Cffs = sdpvar(length(Monoms),1); % Coefficients

% EnerV = dot(Monoms,Cffs); % Lyapunov function


[EnerV,Cffs,Monoms] = polynomial(dX,2,2);
 dVdt = dot(jacobian(EnerV,dX),ddXdt); % Time derivative of Lyapunov function
%% Constraints and Solve

optvar = Cffs; % Let solver tweak the coefficients

EnerVsos = sos(EnerV-epsln*dot(dX,dX));
negdVdtsos = sos(-dVdt-epsln*dot(dX,dX));

constr = [EnerVsos, negdVdtsos, Xtru(1)^2 <= 2.3365^2, Xtru(2)^2 <= 3.263^2, Xtru(3) >= 0, Xtru(3) <= 1.7192]; 
options = sdpsettings('solver', 'mosek');
[sol,v,Q,res] = solvesos(constr, [], options, optvar);

DCMP = sosd(EnerVsos);
DCMPcln = clean(DCMP,1e-6);


sosd(EnerVsos) % Display function and derivative in SOS decomposition
sdisplay(sosd(negdVdtsos))
% sdisplay(value(P))








