clear;clc;
yalmip clear

% For simplicity, let's say that x(1),x(2) correspond to P1,P2, and
% x(3),x(4) correspond to q1,q2

x = sdpvar(4,1);

% Create a small epsilon to ensure definiteness!
eps = 0.00001;

% Define the mu variables, the constraint that is put on them, and the
% objective function
mu = sdpvar(2,1);
mupos = [mu(1),mu(2)];
objfunc = mu(1)+mu(2);


% Here are all the derivatives so we can define the semi-definite
% constraints of the energy and its time derivative.

syst = [-(1+mu(1)) 1.5 0.5 0.2;
    1 -(2+mu(2)) 1 1; 
    0.5 0.5 -3 -1; 
    0.5 0.5 1 -4];
partials = syst*x;

%dP1dt = -(1+mu(1))*x(1)+1.5*x(2)+0.5*x(3)+0.2*x(4);
%dP2dt = x(1)-(2+mu(2))*x(2)+x(3)+x(4);
%dq1dt = 0.5*x(1)+0.5*x(2)-3*x(3)-x(4);
%dq2dt = 0.5*x(1)+0.5*x(2)+x(3)-4*x(4);
%partials = [dP1dt;dP2dt;dq1dt;dq2dt];


Terms = [monolist(x,4,4);monolist(x,2,2)];
Coeffs = sdpvar(length(Terms),1);

% Define the tweakable variables to be the coefficients!
optvar = Coeffs;


% Define the energy function and its positive constraint
U = Coeffs'*Terms;
Usos = U-eps*(x'*x);


% Define the time-derivative of the energy function and its negative
% constraint.
dUdt = jacobian(U,x)*partials;
negdUdtsos = -dUdt-eps*(x'*x);

% When I add in the constraint that sos(negdUdtsos), I receive the error 
% "mosek does not support quadratic semidefinite constraints". I'll try to
% figure out why only the negative derivative is an issue. Removing it
% allows MOSEK to carry out the optimization easily.

cnstr = [sos(Usos), sos(negdUdtsos), mupos(1) >= 0, mupos(2) >= 0]

options = sdpsettings('solver', 'mosek')%, 'sos.congruence', 1, 'sos.reuse', 1,'sos.newton', 1,'sos.numblkdg',1e-4);

%[sol,v,M,res] = 
solvesos(cnstr,objfunc,options,optvar);

% sdisplay(value(mu))

