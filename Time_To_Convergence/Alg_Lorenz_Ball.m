% This script is meant to calculate the minimal nudging coefficients for
% the Lorenz system when we know that the system will at some point be
% within an absorbing ball (because the Lorenz '63 system is dissipative).
clear; clc; yalmip clear;

negdefconst = 1e-26; % Arbitrary constant to ensure negative-definiteness (may not be necessary)


% Define related constants
sig = 10;
beta = 8/3;
rho = 28;

whichcase = [[1;1;1],[1;1;0],[1;0;1],[0;1;1],[1;0;0],[0;1;0]];
MuVals = zeros(3,size(whichcase,2));
% Here's all the YALMIP stuff
for II = 1:size(whichcase,2)
    yalmip clear;

    % Define the decision variables
    mu = sdpvar(3,1);
    X = sdpvar(3,1);
    Xtrue = sdpvar(3,1);

    optvar = mu; % Set mu as the tuned variable

    sym_matrix = [-sig (sig+rho-Xtrue(3))/2 Xtrue(2)/2;
                  (sig+rho-Xtrue(3))/2 -1 0;
                  Xtrue(2)/2 0 -beta]-diag(mu)*diag(whichcase(:,II));
    dUdt = X'*sym_matrix*X
    
    neg_dUdt_sos = -dUdt-negdefconst*dot(X,X); % Negative definiteness constraint
    
    obj_func = dot(mu,mu);
    
    % Introduce the constraints
    constr = [sos(neg_dUdt_sos), mu(:) >= 0, Xtrue(1)^2 <= 2.3365^2, Xtrue(2)^2 <= 3.263^2, Xtrue(3) >= 0, Xtrue(3) <= 1.7192];
    
    options = sdpsettings('solver', 'mosek','verbose',0);
    
    [sol,v,M,res] = solvesos(constr,obj_func,options,optvar);

    
    
    MuVals(:,II) = value(mu);
    
    
end


format long g;
value(MuVals)
format default;