% Runs several simulations of pendulums with time-dependent forcing. Plots
% surfaces showing how the time to convergence varies with nudging
% coefficient and timestep.

%%
clear; clc;
epsilon = 1e-16; % Desired precision
mass = 1.5; % Define important constants
grav = 9.81;
dco = 1.2; 
strl = 2; % Meters, of course.
rtorq = 2.5;
angvel = 0.2;

mumin = mass*strl*(1+grav/strl)^2/(4*dco);

dt = 1e-5; % Initial timestep

% Arbitrary initial conditions
theta01 = 1.3*pi; %.3*pi;
sigma01 = -2.1; %1.2; % Sigma represents the angular velocity
x01 = [theta01;sigma01]; % Put in matrix form

theta02 = 0; %-.4*pi;
sigma02 = 0; %.5;
x02 = [theta02;sigma02]; % Put in matrix form

dtheta0 = theta02-theta01; % Initial difference values for angle and angular velocity
dsigma0 = sigma02-sigma01;

MU = diag([mumin,0]); % Diagonal matrix for MU.
Opsize = size(MU,1);

Oper = @(x,~,t) [x(2);-grav/strl*sin(x(1))-dco/(mass*strl)*x(2)+rtorq/(mass*strl^2)*cos(angvel*t)];
%{
%% Compute Linearized TTC Bound
% Linearized matrix to compute an estimated time to convergence
lin_matrix = [-MU(1,1) 1; 
             grav/strl -dco/(mass*strl)];

[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[dtheta0;dsigma0];
Discrete_Time = log(epsilon/(norm(P_linear)*norm(z0)))/log(norm(eye(Opsize)+D_linear*dt))*dt;
% The nested diag() functions are to prevent the exp() function from
% producing zeros where there should be none. 

%% Compute Analytic TTC Bound

Analytic_Time = log(epsilon/(norm(z0)*norm(P_linear)))/log(norm(expm(D_linear)));
display(Analytic_Time)

max_time = max([Analytic_Time,])

x01 = [theta01;sigma01];
x02 = [theta02;sigma02];

pend1 = FrdEulCntrlTD(x01,Oper,max_time,dt);
pend2 = FrdEulCntrlTD(x02,Oper,max_time,dt,pend1,MU);

pendiff = pend2-pend1;
%%
TTC = TTCfinder(pendiff,dt);
TTC_max = max(TTC);


%% Plotting!
timeline = (0:size(pend1,2)-1)*dt;

%}

%% Figure-Making Code
TIME = 250;
NumTS = 9; % Number of timesteps
rcoff = (1e-1/1e-5)^(1/(NumTS-1)); % Geometrically distributed between 1e-5 and 0.1
Timestep = transpose((1e-5)*rcoff.^(0:NumTS-1)); % Range of timesteps to compute over

EE = length(Timestep);
HH = 50; % However many mu values
barr = 0; % Arbitrary bottom limit
MUmax = mumin;
MUstep = (barr:(MUmax-barr)/(HH-1):MUmax); % Linear spacing between barr and mumin.

[FE_Atime FE_Dtime FE_Rtime Convgs CFL] = deal(zeros(EE,HH));


dtheta0 = theta02-theta01; % All initial conditions we'll need
dsigma0 = sigma02-sigma01;
dx0 = [dtheta0;dsigma0];
x01 = [theta01;sigma01];
x02 = [theta02;sigma02];

for II = 1:EE
    for JJ = 1:HH
        if MUstep(JJ)<2/Timestep(II) % Don't violate CFL condition
            dt = Timestep(II);
            MU = diag([MUstep(JJ) 0]); % Set timestep and MU matrix
            
            data1 = FrdEulCntrlTD(x01,Oper,TIME,Timestep(II),[],[]); % Run reference model
            pred1 = FrdEulCntrlTD(x02,Oper,TIME,Timestep(II),data1,MU); % Run DA model
            pendiff = pred1-data1; % Difference matrix to input to find TTC
        
            FE_Rtime(II,JJ) = max(TTCfinder(pendiff,dt)); % When all variables agree (if they ever do)

            if isnan(pendiff(1,1)) % If pendiff is NaN, because one of the two was too big, no convergence.
                FE_Rtime(II,JJ) = NaN;
            end
    
            if ~(FE_Rtime(II,JJ)<size(pendiff,2)*Timestep(II)) % If it doesn't converge, add in NaN.
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
ylabel('Timestep')
xlabel('MU Value')
zlabel('Numerical TTC')
set(gca,'YScale','log')
