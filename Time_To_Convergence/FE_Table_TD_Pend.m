% Runs several simulations of pendulums with time-dependent forcing. Plots
% surfaces showing how the time to convergence varies with nudging
% coefficient and timestep.
clear; clc;
epsilon = 1e-16; % Desired precision
mass = 1.5; % Define important constants
grav = 9.81;
dco = 1.2; 
strl = 2; % Meters, of course.
rtorq = 250;
angvel = 0.2;

mumin = mass*strl*(1+grav/strl)^2/(4*dco); 
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

MU = diag([mumin,0]); % Diagonal matrix for MU.

operator = @(x,~,t) [x(2);-grav/strl*sin(x(1))-dco/(mass*strl)*x(2)+rtorq/(mass*strl^2)*cos(angvel*t)];

% Linearized matrix to compute an estimated time to convergence
lin_matrix = [-MU(1,1) 1; 
             gravity/strlenth -dampco/(mass*strlenth)];

[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[sd_theta0;sd_sigma0];
% The nested diag() functions are to prevent the exp() function from
% producing zeros where there should be none. 

an_convg_time_lin = log(epsilon/(norm(z0)*norm(P_linear)))/log(norm(expm(D_linear)));
display(an_convg_time_lin)

pend1 = [theta01;sigma01];
pend2 = [theta02;sigma02];

%pend1 = FrdEulCntrlTD(pend1,operator,)

absdiff = abs(pend2-pend1);

TTC = zeros(size(MU,1),1);
for row  = 1: size(TTC,1)
    TTC(row,1) = find(absdiff(row,:),1,"last");
end
TTC_max = max(TTC);

