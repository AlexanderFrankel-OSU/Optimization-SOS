% This script is meant to run a simulation of a forced damped pendulum, and
% see that we can never get synchronization between different initial
% conditions.
clear; clc;

%%
m = 10;
g = 9.81;
b = 1;
L = 1;

dt = 1e-5;
eps = 1e-16;

theta01 = .3*pi;
sigma01 = 1.2;
x01 = [theta01;sigma01];

theta02 = -.4*pi;
sigma02 = .5;
x02 = [theta02;sigma02];

period = pi;
freq_ang = 2*pi/period;
force = 10;

operator = @(x) [x(2);-g/L*sin(x(1))-b/(m*L)*x(2)-force/(m*L^2)*cos(freq_ang*x(1))];
operator2 = @(x,y) [x(2);-g/L*sin(x(1))-b/(m*L)*x(2)-force/(m*L^2)*cos(freq_ang*y(1))];

% Note, you've got to incorporate in the 'true angle' for the forcing to
% work. 

%time = 300;
%pend1 = Forward_Euler(x01,operator,time,dt,[],[]);
%pend2 = Forward_Euler(x02,operator,time,dt,[],[]);


%{
%% Plotting!
absdifftheta12 = abs(pend2(1,:)-pend1(1,:));
plot([0:size(pend1,2)-1]*dt,absdifftheta12) 
% Check that there's no synchronization (there abolutely should not be, 
% since we started with different initial conditions, and the pendulum is
% forced).
%}

%% Data Assimilation Time! We'll use pend1 as the "true system".
MU = diag([300,0]);
Opsize = size(MU,1);


EE = 2;
HH = 2;
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
        MU = diag([300/HH*j,0]);
        MUstep(i,j) = MU(1,1);
        Timestep(i,j) = dt*i;
    
        %% Calculating Analytic Time to Convergence
        [P_dpend,D_dpend] = eig([-MU(1,1) 1; -g/L -b/(m*L)]);
        Analytic_Time  = log(eps/(norm([1 0]*P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(diag(exp(diag(D_dpend)))));
        FE_Atime(i,j) = Analytic_Time;
    
    
        %% Calculating Discrete Time to Convergence
        Discrete_Time = log(eps/(norm([1 0]*P_dpend)*norm(P_dpend^(-1)*sd_x0)))/log(norm(eye(Opsize)+D_dpend*Timestep(i,j)))*Timestep(i,j);
        FE_Dtime(i,j) = Discrete_Time;
    
    
        %% Calculating Actual Time to Convergence According to MATLAB
        max_time = max([Analytic_Time,Discrete_Time])+20;

        x01 = [theta01;sigma01];
        x02 = [theta02;sigma02];

        pend1 = Forward_Euler(x01,operator,max_time,Timestep(i,j),[],[]);
        prediction = Forward_Euler(x02,operator2,max_time,Timestep(i,j),pend1,MU);
        absdifftheta12 = abs(pend1(1,:)-prediction(1,:));
    
    
        for Num = 1:size(absdifftheta12,2) 
            if absdifftheta12(1,Num)<eps
                break;
            end
        end
        Num = Num+1;
        Real_Time = Num*dt;
        FE_Rtime(i,j) = Real_Time;
    
    end
end

%% PLOT
surf(Timestep,MUstep,transpose(FE_Rtime))
xlabel('Timestep')
ylabel('MUstep')
zlabel('FE Numerical TTC')