clear;clc;
yalmip clear

% Create a single sdpvar with three elements to call
x = sdpvar(3,1);

% Create sdpvars that can represent the terms and coeffcients of the
% polynomial

Terms = monolist([x(1),x(2),x(3)],2);
Coeffs = sdpvar(length(Terms),1);
% Define the optvar that MOSEK will be tuning (it shouldn't be tuning the
% terms, Xander)
optvar = Coeffs;

% Define V by the sum of the product of the terms and coefficients
V = sum(Coeffs.*Terms);

% Create a small epsilon to turn semidefinite constraints in to definite
eps = 0.00001;
sig = 10; r = 0.5; b = 8/3;

% Define the derivatives with respect to time, given by the Lorenz system
dxdt = sig*(x(2)-x(1));
dydt = x(1)*(r-x(3))-x(2);
dzdt = x(1)*x(2)-b*x(3);

% Define the derivative of V using the jacobian function
dVdt = sum(jacobian(V,x(1))*dxdt+jacobian(V,x(2))*dydt+jacobian(V,x(3))*dzdt);

% Define specific constrants
% Constraint 1: the coefficients of V should make it sum of squares
Vsos = sum(Coeffs.*Terms)-eps*(x(1)^2+x(2)^2+x(3)^2);
% Constraint 2: If we want the derivative of V to be strictly negative, 
% we hold that the negative of it should be a sum of squares, i.e., 
% positive definite.
negdVdt_sos = -dVdt-eps*(x(1)^2+x(2)^2+x(3)^2); 

F_cnstr = [sos(Vsos),  sos(negdVdt_sos)];

% Here we set the options for YALMIP to use MOSEK, and what valiables
% should receive the outputs
options = sdpsettings('solver', 'mosek');
[sol,v,Q,res] = solvesos(F_cnstr, [], options, optvar);

% Creating a coefficient vector that uses the values of the Coeffs list of
% terms
coeff_vec = cell(length(Coeffs),1);
for i = 1:length(Coeffs)
    coeff_vec{i} = value(Coeffs(i));
end

% Doing the same for the polynomial itself, 'plnm' 
plnm = 0;
for i = 1:length(Terms)
plnm = plnm + coeff_vec{i}*Terms(i);
end

disp('Sum of squares Lyapunov function:')
sdisplay(plnm)
disp('Derivative of Lyapunov function')
sdisplay(.5*(jacobian(plnm,x(1))*dxdt+jacobian(plnm,x(2))*dydt+jacobian(plnm,x(3))*dzdt))
disp('These are the terms and their corresponding coefficients')
sdisplay([value(Coeffs),Terms])
