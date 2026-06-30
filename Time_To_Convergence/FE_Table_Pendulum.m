clear; clc;
% Here, sd is short for lowercase delta.
% This script may not be used, as it isn't a great example of data
% assimilation, since the system decays to a single state anyway.
m = 10; % Useful constants (may be interesting to run through many possible values of b)
g = 9.81;
b = 1;
dt = 1e-5;
eps = 1e-16;
MU = diag([1,0]);

theta0t = .9*pi; % Defining all initial conditions
theta0 = -.8*pi;

theta0t = mod(theta0t,2*pi);
theta0 = mod(theta0,2*pi);

sd_theta0 = theta0-theta0t;
sigma0t = 1.5;
sigma0 = -1;
sd_sigma0 = sigma0-sigma0t;
x0t = [theta0t;sigma0t];
x0 = [theta0;sigma0];

%% Computing Linearized Analytic Time to Convergence
lin_matrix = [-MU(1,1) 1; % Linearized matrix to compute an estimated time to convergence
             -g -b/m];
[P_linear,D_linear] = eig(lin_matrix);
z0 = (P_linear^-1)*[sd_theta0;sd_sigma0];

an_convg_time_lin = log(eps/(norm(z0)*norm(P_linear)))/log(norm(diag(exp(diag(D_linear)))));
display(an_convg_time_lin)



operator = @(x) [x(2);-g*sin(mod(x(1),2*pi))-b/m*x(2)];
FE_data_pend = Forward_Euler(x0t,operator,4*an_convg_time_lin,dt,[],[]);
FE_pred_pend = Forward_Euler(x0 ,operator,4*an_convg_time_lin,dt,FE_data_pend,MU);
FE_pend_absdiff = abs(FE_pred_pend-FE_data_pend);


for num = 1:size(FE_pend_absdiff,2)
    if FE_pend_absdiff(1,num)<eps
        break;
    end
end
num = num+1;
FE_disc_pend_TTC = num*dt;
display('Time to convergence according to MATLAB:')
display(FE_disc_pend_TTC)

%% Run the Itterations
