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
clear; clc;

%%
epsilon = 1e-16; % Desired precision

mass = 1.5; % Define important constants
gravity = 9.81;
dampco = 1.2; 
strlenth = 2; % Meters, of course.
torque = 250; % Arbitrary torque value, may be a bit like attaching a small rocket engine to the pendulum
% If the applied torque is less than mgL then the system has an equilibrium

% The minimum mu value can be calculated analytically because this is a
% two-variable system
mumin = mass*strlenth*(1+gravity/strlenth)^2/(4*dampco); 

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



MU = diag([10,0]); % Diagonal matrix for MU.

operator = @(x) [x(2);-gravity/strlenth*sin(x(1))-dampco/(mass*strlenth)*x(2)+torque/(mass*strlenth^2)]; % Define the operator according to the equation, where x(2) is angular velocity
%{
SGN = @(z) (1)*(z>=0)+(-1)*(z<0);
direction = @(x,y) (SGN(sin(y(1)-x(1))));
angdist = @(x,y) (min(mod(x(1)-y(1),2*pi),2*pi-mod(x(1)-y(1),2*pi)));
MUprtr = @(x,y,unnec3) [direction(x(1),y(1))*mumin*angdist(x(1),y(1));0];
%% Calculating Bounds on (Linearized) Analytic and Discrete Time to Convergence

% Linearized matrix to compute an estimated time to convergence
lin_matrix = [-MU(1,1) 1; 
             gravity/strlenth -dampco/(mass*strlenth)];

[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[sd_theta0;sd_sigma0];
% The nested diag() functions are to prevent the exp() function from
% producing zeros where there should be none. 

an_convg_time_lin = abs(log(epsilon/(norm(z0)*norm(P_linear)))/log(norm(expm(D_linear))));
display(an_convg_time_lin)
Opsize = size(MU,1);


pend1 = FrdEulCntrlTD(x01,operator,an_convg_time_lin,dt); % Create the "data"
pend2 = FrdEulCntrlTD(x02,operator,an_convg_time_lin,dt); % Run the model with different initial conditions and apply DA

pendiff = pend2-pend1;

TTC = zeros(Opsize,1);
            if anynan(pendiff(:,end)) % Calculate the exact TTC, by whatever means necessary.
                 notNAN = zeros(Opsize,1);
                 startNAN = size(pendiff,2);
                 for row = 1:Opsize
                    notNAN(row,1) = find(~isnan(pendiff(row,1:startNAN)),1,"last");
                    startNAN = notNAN(row,1);
                 end
                 notNANmin = min(notNAN);
                 for row = 1:Opsize
                    TTC(row,1) = (find(pendiff(row,1:notNANmin),1,"last")+1)*dt;
                 end
            else
                for row = 1:Opsize
                    TTC(row,1) = (find(pendiff(row,:),1,"last")+1)*dt;
                end
            end
            
    TTC_max = max(TTC)

%% Plotting!
timeline = (0:size(pend1,2)-1)*dt; % Useful for creating plots.
absdiff12 = abs(pend2-pend1); 
semilogy(timeline,absdiff12);
% Check that there's no synchronization (there should not be, 
% since we started with different initial conditions, and the pendulum is
% forced).

%}

%% Table-Making Code

%  We'll use pend1 as the "true system" that pend2 is driven to chase by
%  the nudging coefficient.

NumTS = 11;
rcoff = (1e+5)^(1/(NumTS-1));
Timestep = (1e-5)*rcoff.^(0:NumTS-1); % Range of timesteps to compute over

EE = size(Timestep,2); 
HH = 50; % However many MU values
barr = 0.01; % Arbitrary bottom limit
MUmax = mumin;
MUstep = (barr:(MUmax-barr)/(HH-1):MUmax); % Linear spacing between barr and mumin.

MU = diag([mumin,0]);
Opsize = size(MU,1);

FE_Atime = zeros(EE,HH);
FE_Dtime = zeros(EE,HH);
FE_Rtime = zeros(EE,HH);
Convgs = zeros(EE,HH);
NIT = ones(EE,HH);

sd_theta0 = theta02-theta01;
sd_sigma0 = sigma02-sigma01;
sd_x0 = [sd_theta0;sd_sigma0];

x01 = [theta01;sigma01];
x02 = [theta02;sigma02];

for II = 1:EE
    for JJ = 1:HH
        
        dt = Timestep(1,II);
        MU = diag([MUstep(JJ) 0]);
    
        %% Calculating Analytic Time to Convergence
        [P_dpend,D_dpend] = eig([-MU(1,1) 1; gravity/strlenth -dampco/(mass*strlenth)]);
        Analytic_Time  = abs(log(epsilon/(norm(P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(expm(D_dpend))));
        FE_Atime(II,JJ) = Analytic_Time;
        
    
        %% Calculating Discrete Time to Convergence
        Discrete_Time = abs(log(epsilon/(norm(P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(eye(Opsize)+D_dpend*Timestep(1,II)))*Timestep(1,II));
        FE_Dtime(II,JJ) = Discrete_Time;
    
  
        %% Calculating Actual Time to Convergence According to MATLAB
        while or(Convgs(II,JJ)==0,NIT(II,JJ)==1)
        max_time = NIT(II,JJ)*max([Analytic_Time,Discrete_Time]);
        if NIT(II,JJ) == 2
            max_time = max([max_time,600]);
        end

        data1 = Forward_Euler(x01,operator,max_time,Timestep(1,II),[],[]);
        pred1 = Forward_Euler(x02,operator,max_time,Timestep(1,II),data1,MU);
        pendiff = data1-pred1;
        
            TTC = zeros(Opsize,1);
            if anynan(pendiff(:,end))
                 notNAN = zeros(Opsize,1);
                 startNAN = size(pendiff,2);
                 for row = 1:Opsize
                    notNAN(row,1) = find(~isnan(pendiff(row,1:startNAN)),1,"last");
                    startNAN = notNAN(row,1);
                 end
                 notNANmin = min(notNAN);
                 for row = 1:Opsize
                    TTC(row,1) = (find(pendiff(row,1:notNANmin),1,"last")+1)*dt;
                 end
            else
                for row = 1:Opsize
                    TTC(row,1) = (find(pendiff(row,:),1,"last")+1)*dt;
                end
            end
        FE_Rtime(II,JJ) = max(TTC);
        
        if NIT(II,JJ) == 2
            break;
        end
        if FE_Rtime(II,JJ)<max_time
            Convgs(II,JJ) = 1;
        else
            NIT(II,JJ) = 2;
            Convgs(II,JJ) = 0;
        end
    
        end
    end
end



%% Plots (Feel free to edit these)

surf(Timestep,MUstep,transpose(FE_Rtime))
xlabel('Timestep')
ylabel('MUstep')
zlabel('FE Numerical TTC')


%{
surf(log10(Timestep),MUstep,transpose(FE_Dtime))
xlabel('log10(Timestep)')
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