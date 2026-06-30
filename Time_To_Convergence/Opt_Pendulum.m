% This script is meant to calculate the minimal nudging coefficient for a
% damped pendulum, and then use this to calculate the time to convergence
% for the two numerical schemes we're using: Forward Euler and
% Adams-Bashforth.
clear; clc; yalmip clear

%% Run YALMIP to Find Minimal MU
m = 10; % Useful constants (may be interesting to run through many possible values of b)
g = 9.81;
b = 1;

z = sdpvar(2,1); % Define the two variables
mu = sdpvar(1);
optvar = mu; % Define the value the solver tries to tweak
eps = 1e-9; % Add small positive constant to ensure negative-definiteness
Q = [-mu 0.5*(1+g); % Write the symmetric matrix in terms of the variables and constants
    0.5*(1+g) -b/m];
dUdt = z'*Q*z

neg_dUdt_sos = -dUdt-eps*dot(z,z); % Negative definiteness constraint

obj_func = mu^2;

constr = sos(neg_dUdt_sos);

options = sdpsettings('solver', 'mosek');

[sol,v,M,res] = solvesos(constr,obj_func,options,optvar);

sdisplay(value(mu))
