% This script is meant to make a table of values from calculating many
% simulations of a forced damped pendulum with data assimilation. The
% forcing term is a constant (it's really a rtorq).
% The system is derived from the rtorq-balance equation:
% (mL^2)0'' = -mgLsin(0)-bL0'+T
% Here "0" stands for angle from straight down, b is the drag coefficient,
% and T is the constant forcing (torquing) term in the counterclockwise direction.

%% Corresponding Figure/Table/Values in Paper
% This script produces figures and values corresponding to the section in
% the paper regarding the forced damped pendulum.

%%
clear; clc;

epsilon = 1e-16; % Desired precision

mass = 1.5; % Define important constants
grav = 9.81;
dco = 1.2; 
strl = 2; % Meters, of course.
rtorq = 250; % Arbitrary rtorq value, may be a bit like attaching a small rocket engine to the pendulum
% If the applied rtorq is less than mgL then the system has an equilibrium

% The minimum mu value can be calculated analytically because this is a
% two-variable system
mumin = mass*strl*(1+grav/strl)^2/(4*dco); 

dt = 1e-5; % Initial timestep


% Arbitrary initial conditions
theta01 = .3*pi;
sigma01 = 1.2; % Sigma represents the angular velocity
x01 = [theta01;sigma01]; % Put in matrix form

theta02 = 0; %-.4*pi;
sigma02 = 0; %.5;
x02 = [theta02;sigma02]; % Put in matrix form


dtheta0 = theta02-theta01; % Initial difference values for angle and angular velocity
dsigma0 = sigma02-sigma01;

MU = diag([5,0]); % Diagonal matrix for MU.


operator = @(x) [x(2);-grav/strl*sin(x(1))-dco/(mass*strl)*x(2)+rtorq/(mass*strl^2)]; % Define the operator according to the equation, where x(2) is angular velocity
%{
%% Calculating Bounds on (Linearized) Analytic TTC

% Linearized matrix to compute an estimated time to convergence
lin_matrix = [-MU(1,1) 1; 
             grav/strl -dco/(mass*strl)];

[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[dtheta0;dsigma0];

an_convg_time_lin = abs(log(epsilon/(norm(z0)*norm(P_linear)))/log(norm(expm(D_linear))));
display(an_convg_time_lin)
Opsize = size(MU,1);


%% Calculating Bound on (Linearized) Discretized TTC
Discrete_Time = abs(log(epsilon/(norm(P_linear)*norm(z0)))/log(norm(eye(Opsize)+D_linear*dt))*dt);
display(Discrete_Time)


%%
max_time = max([an_convg_time_lin,Discrete_Time]);
pend1 = FrdEulCntrlTD(x01,operator,max_time,dt); % Create the "data"
pend2 = FrdEulCntrlTD(x02,operator,max_time,dt,pend1,MU); % Run the model with different initial conditions and apply DA

pendiff = pend2-pend1;
TTC = TTCfinder(pendiff,dt);
TTC_max = max(TTC)

%% Plotting!
timeline = (0:size(pend1,2)-1)*dt; % Useful for creating plots.
absdiff12 = abs(pend2-pend1); 
semilogy(timeline,absdiff12);
% plot(timeline,pend1)
% Check that there's no synchronization (there should not be, 
% since we started with different initial conditions, and the pendulum is
% forced).

%}

%% Table-Making Code

%  We'll use pend1 as the "true system" that pend2 is driven to chase by
%  the nudging coefficient.

NumTS = 9;
rcoff = (1e-1/1e-5)^(1/(NumTS-1));
Timestep = transpose((1e-5)*rcoff.^(0:NumTS-1)); % Range of timesteps to compute over

EE = size(Timestep,1); 
HH = 50; % However many MU values
barr = 0; % Arbitrary bottom limit
MUmax = mumin;
MUstep = (barr:(MUmax-barr)/(HH-1):MUmax); % Linear spacing between barr and mumin.

MU = diag([mumin,0]);
Opsize = size(MU,1);

[FE_Atime FE_Dtime FE_Rtime Convgs] = deal(zeros(EE,HH));
CLF = ones(EE,HH);

dtheta0 = theta02-theta01;
dsigma0 = sigma02-sigma01;
dx0 = [dtheta0;dsigma0];

x01 = [theta01;sigma01];
x02 = [theta02;sigma02];

for II = 1:EE
    for JJ = 1:HH
        
        dt = Timestep(II);
        MU = diag([MUstep(JJ) 0]);
    
        %% Calculating Linearized Analytic Time to Convergence
        [P_dpend,D_dpend] = eig([-MU(1,1) 1; grav/strl -dco/(mass*strl)]);
        Analytic_Time  = log(epsilon/(norm(P_dpend)*norm(P_dpend^(-1)*dx0)))/log(norm(expm(D_dpend)));
        FE_Atime(II,JJ) = Analytic_Time;
        
    
        %% Calculating Linearized Discrete Time to Convergence
        Discrete_Time = log(epsilon/(norm(P_dpend)*norm(P_dpend^(-1)*dx0)))/log(norm(eye(Opsize)+D_dpend*Timestep(II)))*Timestep(II);
        FE_Dtime(II,JJ) = Discrete_Time;
    
        %% Calculating Actual Time to Convergence According to MATLAB
        
        TIME = 250;
        if MUstep(JJ)<2/Timestep(II)
            data1 = FrdEulCntrlTD(x01,operator,TIME,Timestep(II));
            pred1 = FrdEulCntrlTD(x02,operator,TIME,Timestep(II),data1,MU);
            pendiff = pred1-data1;
    
            FE_Rtime(II,JJ) = max(TTCfinder(pendiff,dt)); % Exact convergence happens here
            
            if isnan(pendiff(1,1)) % If pendiff is NaN, because one of the two was too big, no convergence.
                FE_Rtime(II,JJ) = NaN;
            end
    
            if ~(FE_Rtime(II,JJ) < TIME) % If it doesn't converge, add in NaN.
                FE_Rtime(II,JJ) = NaN;
            else % If actual convergence happens
                Convgs(II,JJ) = 1;
            end
        else
            FE_Rtime(II,JJ) = NaN;
            CFL(II,JJ) = 0;
        end

    end
end


%% Plot!

clf;
surf(MUstep,Timestep,FE_Rtime)
ylabel('Step Size')
xlabel('Nudging Coefficient')
zlabel('Numerical TTC')
set(gca,'YScale','log')





%% This is not to be run. It just describes how you could find the same value of mu using SOS Optimization.

% clear; clc; yalmip clear
% z = sdpvar(2,1); % Define the two variables
% mu = sdpvar(1);
%
% optvar = mu; % Define the value the solver tries to tweak
%
% negdefconst = 1e-9; % Add small positive constant to ensure negative-definiteness
%
% Q = [-mu 0.5*(1+g/L); % Write the symmetric matrix in terms of the variables and constants
%     0.5*(1+g/L) -b/(m*L)];
% dUdt = z'*Q*z % Define derivative of energy Lyapunov function
% 
% neg_dUdt_sos = -dUdt-negdefconst*dot(z,z); % Negative definiteness constraint
% 
% obj_func = mu^2; % Minimizing mu
% 
% constr = [sos(neg_dUdt_sos),mu >= 0];
% 
% options = sdpsettings('solver', 'mosek');
% 
% [sol,v,M,res] = solvesos(constr,obj_func,options,optvar);
% 
% sdisplay(value(mu))