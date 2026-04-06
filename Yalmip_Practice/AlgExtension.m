clear;clc;
yalmip clear

% For simplicity, let's say that x(1),x(2) correspond to P1,P2, and
% x(3),x(4) correspond to q1,q2

x = sdpvar(4,1);

% Create a small epslon to ensure definiteness!
eps = 0.00001;

% Define the mu variables, the constraint that is put on them, and the
% objective function
mu = sdpvar(2,1);
mupos = [mu(1)-eps,mu(2)-eps];
objfunc = mu(1)+mu(2);




% Here are all the derivatives so we can define the semi-definite
% constraints of the energy and its time derivative.

dP1dt = -(1+mu(1))*x(1)+1.5*x(2)+0.5*x(3)+0.2*x(4);
dP2dt = x(1)-(2+mu(2))*x(2)+x(3)+x(4);
dq1dt = 0.5*x(1)+0.5*x(2)-3*x(3)-x(4);
dq2dt = 0.5*x(1)+0.5*x(2)+x(3)-4*x(4);
Partials = [dP1dt,dP2dt,dq1dt,dq2dt]';

Terms = monolist(x,4,2);
Coeffs = sdpvar(length(Terms),1);

% Define the tweakable variables to be the coefficients!
optvar = [mu;Coeffs];


% Define the energy function and its positive constraint
U = Coeffs'*Terms;
Usos = U-eps*(x'*x);


% Define the time-derivative of the energy function and its negative
% constraint.
dUdt = sum(jacobian(U,x)*Partials);
negdUdtsos = -dUdt-eps*(x'*x);


Constr = [sos(Usos), sos(negdUdtsos), mupos(1) >= 0, mupos(2) >= 0]

options = sdpsettings('solver', 'mosek','sos.congruence',1, 'sos.reuse', 0);

[sol,v,M,res] = solvesos(Constr,objfunc,options,optvar);

sdisplay(value(mu))

