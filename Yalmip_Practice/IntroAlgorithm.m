clear;clc;
yalmip clear

% Define sdpvars x and mu, and set mu to be the optimization variable
% (we're not tweaking the p_s or q_s [which are the four entries of x]).
x = sdpvar(4,1);
mu = sdpvar(2,1);
optvar = mu;
% Create a small epsilon to ensure definiteness!
eps = 1e-9;

% Define Q as a symmetric matrix (by turning a triangular sum of quadratic
% terms into the form z'*Q*z from the system).
KUE = [2 4.25 0.5 0.2;
            1 3 1 1; 
            0.5 0.5 -3 -1; 
            0.5 0.5 1 -4];
Q = (KUE+KUE')/2-diag([mu(1),mu(2),0,0]);


% Define the energy term (as far as I know, the 0.5 term isn't significant
% for the optimization; though it caused the solver to give slightly
% different values for the mu(2) parameter

dUdt = x'*Q*x;
%Subtract a small fraction of the sum of squares to ensure definiteness
negdUdtsos = -dUdt-eps*(dot(x,x));

% Define other constraints, such as that
objfunc = dot(mu,mu);

% Define the constraints on the optimization. The mu_s need to be strictly
% positive, and Q should be negative semi-definite (which is equivalent to
% the expression for dUdt being negative semi-definite, which MOSEK can
% handle).
constr = [sos(negdUdtsos), mu(:) >= 0];

options = sdpsettings('solver', 'mosek');

[sol,v,M,res] = solvesos(constr,objfunc,options,optvar);

sdisplay(value(mu))



