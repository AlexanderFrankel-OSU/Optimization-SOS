% This script is meant to find Lyapunov Functions for particular systems of
% differential equations.

% Currently does not work.
clear; clc; yalmip clear;
grav = 9.81;strl = 2;rtorq = 250;dco = 1.2;mass = 1.5; % Constants

MOmin = mass*strl*(1+grav/strl)^2/(4*dco); % Analytically calculated from algorithm
mu = 20; % Just for this pendulum, we'll change this
MU = diag([mu 0 0 0]);
muvec = diag(MU); % Useful for defining the difference derivatives dsXdt.

epsln = 1e-8;
LinOp = [-1 2 0.5 0.2;
             1 -2 1 1;
             0.5 0.5 -3 -1;
             0.5 0.5 1 -4];
Oper = @(x) (LinOp*x)
%[x(2)-mu*x(1);-grav/strl*sin(x(1))-dco/(mass*strl)*x(2)+rtorq/(mass*strl^2)];

numvars = 4;

Xtru = sdpvar(numvars,1); % The true system
X = sdpvar(numvars,1); % A model of the system
sX = sdpvar(numvars,1); % Difference variables
dXtrudt = Oper(Xtru);
dXdt = Oper(X)-MU*(sX);
%dsXdt = dXdt-dXtrudt; % Difference operator with DA (if any)

dsXdt = Oper(sX)-MU*sX;

Monoms = monolist(sX,2,2); % Monomial of difference variables
Cffs = sdpvar(size(Monoms,1),1);

optvar = Cffs;

EnerV = dot(Monoms,Cffs);
Vsos = sos(EnerV-epsln*dot(sX,sX));
dVdt = dot(jacobian(EnerV,sX),dXdt);
negdVdtsos = sos(-dVdt-epsln*dot(sX,sX));

constr = [Vsos, negdVdtsos]; 
options = sdpsettings('solver', 'mosek');
[sol,v,Q,res] = solvesos(constr, [], options, optvar);

sdisplay([value(Cffs),Monoms])
sdisplay(dot(value(Cffs),Monoms))

%sdisplay(dot(jacobian(EnerV,X),dXdt)+dot(jacobian(EnerV,Xtru),dXtrudt))









