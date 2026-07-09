% This script is meant to make a table of values from calculating many
% simulations of a forced damped pendulum with data assimilation. The
% forcing term is a constant (it's really a torque).
% The system is derived from the torque-balance equation:
% (mL^2)0'' = -mgLsin(0)-bL0'+T
% Here "0" stands for angle from straight down, b is the drag coefficient,
% and T is the constant forcing (torquing) term in the counterclockwise direction.

%% Corresponding Figure/Table/Values in Paper
% This script produces figures and values corresponding to the section in
% the paper regarding the forced damped pendulum.
display(['This script produces figures and values corresponding to the\ ' ...
         'section in the paper regarding the forced damped pendulum'])


%%
clear; clc;
epsilon = 1e-16; % Desired precision

mass = 1.5; % Define important constants
gravity = 9.81;
dampco = 1.2; 
strlenth = 2; % Meters, of course.

% The minimum mu value can be calculated analytically because this is a
% two-variable system
mumin = mass*strlenth*(1+gravity/strlenth)^2/(4*dampco); 
display(mumin)

dt = 1e-5; % Initial timestep


% Arbitrary initial conditions
theta01 = .3*pi;
sigma01 = 1.2; % Sigma represents the angular velocity
x01 = [theta01;sigma01]; % Put in matrix form

theta02 = -.4*pi;
sigma02 = .5;
x02 = [theta02;sigma02]; % Put in matrix form


sd_theta0 = theta02-theta01; % Initial difference values for angle and angular velocity
sd_sigma0 = sigma02-sigma01;

torque = 250; % Arbitrary torque value, may be a bit like attaching a small rocket engine to the pendulum
% If the applied torque is less than mgL then the system has an equilibrium

MU = diag([mumin,0]); % Diagonal matrix for MU.

operator = @(x) [x(2);-gravity/strlenth*sin(x(1))-dampco/(mass*strlenth)*x(2)+torque/(mass*strlenth^2)]; % Define the operator according to the equation, where x(2) is angular velocity

%% Calculating Bounds on (Linearized) Analytic and Discrete Time to Convergence

% Linearized matrix to compute an estimated time to convergence
lin_matrix = [-MU(1,1) 1; 
             -gravity/strlenth -dampco/(mass*strlenth)];

[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[sd_theta0;sd_sigma0];
% The nested diag() functions are to prevent the exp() function from
% producing zeros where there should be none. 

an_convg_time_lin = log(epsilon/(norm(z0)*norm([1 0]*P_linear)))/log(norm(diag(exp(diag(D_linear)))));
display(an_convg_time_lin)


pend1 = Forward_Euler(x01,operator,an_convg_time_lin,dt,[],[]); % Create the "data"
pend2 = Forward_Euler(x02,operator,an_convg_time_lin,dt,pend1,MU); % Run the model with different initial conditions and apply DA

%% Plotting!
timeline = (0:size(pend1,2)-1)*dt; % Useful for creating plots.
for col = 1:size(pend1,2)
    absdifftheta12 = norm(pend2(:,col)-pend1(:,col)); 
end
semilogy(timeline,absdifftheta12); 
% Check that there's no synchronization (there should not be, 
% since we started with different initial conditions, and the pendulum is
% forced).


%% Table-Making Code
promp = "Continue to Table-Maker? [Y/N] Anything but 'Y' or 'N' is no."
goforth = stringtype(promp,"s")

if isequal(goforth,'Y')

elseif isequal(goforth,'N')
    return;
else
    return;
end

%  We'll use pend1 as the "true system" that pend2 is driven to chase by
%  the nudging coefficient.

MU = diag([mumin,0]);
Opsize = size(MU,1);


EE = 4; % Number of different timesteps to run through, evenly spaces them between dt and EE*dt. Ideally keep this number low.
if EE < 1
    display('You probably want a timestep. At least one...')
    return;
end
HH = 4; % Will cause error if zero.

if HH < 1
    display('You probably want some values for mu. That is the whole point.')
    return;
end 
FE_Atime = zeros(EE,HH);
FE_Dtime = zeros(EE,HH);
FE_Rtime = zeros(EE,HH);
Timestep = zeros(1,EE);
MUstep = zeros(1,HH);

sd_theta0 = theta02-theta01;
sd_sigma0 = sigma02-sigma01;
sd_x0 = [sd_theta0;sd_sigma0];

for i = 1:EE
    for j = 1:HH
        MU = diag([mumin/HH*j,0]);
        MUstep(1,j) = MU(1,1);
        Timestep(1,i) = dt*i;
    
        %% Calculating Analytic Time to Convergence
        [P_dpend,D_dpend] = eig([-MU(1,1) 1; -gravity/strlenth -dampco/(mass*strlenth)]);
        Analytic_Time  = log(epsilon/(norm([1 0]*P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(diag(exp(diag(D_dpend)))));
        FE_Atime(i,j) = Analytic_Time;
        
    
        %% Calculating Discrete Time to Convergence
        Discrete_Time = log(epsilon/(norm([1 0]*P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(eye(Opsize)+D_dpend*Timestep(1,i)))*Timestep(1,i);
        FE_Dtime(i,j) = Discrete_Time;
    
    
        %% Calculating Actual Time to Convergence According to MATLAB
        max_time = max([Analytic_Time,Discrete_Time]);

        x01 = [theta01;sigma01];
        x02 = [theta02;sigma02];

        data1 = Forward_Euler(x01,operator,max_time,Timestep(1,i),[],[]);
        prediction = Forward_Euler(x02,operator,max_time,Timestep(1,i),data1,MU);
        absdifftheta12 = abs(data1(1,:)-prediction(1,:));
    
    
        for Num = 1:size(absdifftheta12,2) 
            if absdifftheta12(1,Num)<epsilon % This is a really terrible way to check for convergence. UPDATE.
                break;
            end
        end
        Num = Num+1;
        Real_Time = Num*dt;
        FE_Rtime(i,j) = Real_Time;
    end
end

%% Plots (Feel free to edit these)
surf(Timestep,MUstep,FE_Rtime)
xlabel('Timestep')
ylabel('MUstep')
zlabel('FE Numerical TTC')

%{
surf(Timestep,MUstep,transpose(FE_Dtime))
xlabel('Timestep')
ylabel('MUstep')
zlabel('FE Discrete TTC')
%}



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