%% Run the algorithm
clear;clc;yalmip clear;
epsln = 1e-8;
Oper = [1.7 4; 2.2 -1.3];

dX = sdpvar(2,1); mu = sdpvar(2,1);
MU = diag(mu);

optvar = mu;
obj_func = dot(mu,mu);

dEnergydt = dX'*((Oper+Oper')/2-MU)*dX; % The algorithm has this Lyapunov function
negdEnergydtSOS = sos(-dEnergydt-epsln*dot(dX,dX)); % Time derivative of Lyapunov function
constr = [negdEnergydtSOS, mu(:) >= 0]; % Energy is already a sum of squares
options = sdpsettings('solver','mosek');

[sol,v,Q,res] = solvesos(constr,obj_func,options,optvar);

display('The negative SOS decomposition for the time derivative of the energy is:')
negdeco = sosd(negdEnergydtSOS);
sdisplay(negdeco)
sdisplay(value(mu))
mumins = value(mu);

%% Find a Lyapunov function for lower values of mu
yalmip clear;

dX = sdpvar(2,1);
mu = [2.4;.9];
MU = diag(mu);

[Enrg,Cffs,Monoms] = polynomial(dX,2,2); % Creat a polynomial for the Lyapunov function
dEnrgdt = dot(jacobian(Enrg,dX),(Oper-MU)*dX);

EnrgSOS = sos(Enrg-epsln*dot(dX,dX)); % SOS constraints
negdEnrgdtSOS = sos(-dEnrgdt-epsln*dot(dX,dX));


constr = [EnrgSOS, negdEnrgdtSOS, mu(:) >= 0];
options = sdpsettings('solver','mosek');
optvar = Cffs;

[sol,v,Q,res] = solvesos(constr,[],options,optvar);

display('Energy SOS decomposition')
sdisplay(sosd(EnrgSOS)) % Display the energy SOS decomposition
display('Time derivative of energy function')
sdisplay(sosd(negdEnrgdtSOS)) % Display the energy's time derivative

% %% Calculating minimal coefficients the bad way
% ITTS = 1000;
% 
% 
% musq = zeros(ITTS);
% 
% for II = 1:ITTS
%     for JJ = 1:ITTS
%         mu = [mumins(1)/ITTS*II;mumins(2)/ITTS*JJ];
%         if and(trace(Oper)-trace(MU)<0,(Oper(1,1)-mu(1))*(Oper(2,2)-mu(2))-Oper(2,1)*Oper(1,2)>0)
%            musq(II,JJ) = dot(mu,mu);
%         else
%            musq(II,JJ) = NaN;
%         end
%     end
% end
% 
% [valmax, inmax] = min(musq(:));
% display(inmax)


