% This script is meant to compare the different approximations for the time
% to convergence of a linear system over various nudging coefficients,
% timesteps, and initial conditions.
clear;clc;

x0t = rand(4,1)*10;%[.1;.1;-.1;.1];
x0 = zeros(4,1); %[-.2;-.2;-.2;.1]; Previously used.


sx0 = x0-x0t;
dt = 1e-5;
epsln = 1e-16;

Oper = [-1 2 0.5 0.2;
             1 -2 1 1;
             0.5 0.5 -3 -1;
             0.5 0.5 1 -4];
x1t = x0t+Oper*x0t;
sz1 = [x1t;x0t];
Opsize = size(Oper,1);

%% Begin long for loop

NumTS = 200;
stepend = 1e-2;
stepstart = 1e-7;
rcoff = (stepend/stepstart)^(1/(NumTS-1));
Timestep = transpose((stepstart)*rcoff.^(0:NumTS-1)); % This defines the spread of the timesteps -- Geometrically from 1e-5 to 1e-1.

EE = size(Timestep,1);
HH = 50;
mumin = 0.64;
mumax = 200; % Otherwise run this for the larger, less-interesting plot.
MUstep = (mumin:(mumax-mumin)/(HH-1):mumax);

[FE_Atime FE_Dtime AB_Dtime] = deal(zeros(EE,HH));


for II = 1:EE
    for JJ = 1:HH % Start at zero so that we see no convergence occurs if mu is zero
            
            dt = Timestep(II);
            MU = diag([MUstep(1,JJ),0,0,0]);
            
            
            %% Computing the Analytic Time to Convergence
            
            [P_FE,D_FE] = eig(Oper-MU); % Calculate the eigenvalues
            FE_Atime(II,JJ) = log(epsln/(norm(P_FE)*norm((P_FE^-1)*sx0)))/log(norm(expm(D_FE))); % Calculate upper bound for analytic time to convergence
   
            %% (Forward Euler) Computing the Discrete Time to Convergence
            FE_Dtime(II,JJ) = dt/log(norm(eye(Opsize)+D_FE*dt))*log(epsln/((norm(P_FE)*norm((P_FE^-1)*sx0))));

            %% (Adams-Bashforth) Computing the Discrete Time to Convergence
            G = (eye(Opsize)+(1.5*Oper-MU)*dt);
            H = -0.5*Oper*dt;
            B = [G H;eye(Opsize) zeros(Opsize)];
            [P_AB D_AB] = eig(B);
            AB_Dtime(II,JJ) = dt/log(norm(D_AB))*log(epsln/(norm([eye(Opsize) zeros(Opsize)]*P_AB)*norm(P_AB^-1*sz1)));
    end
end

%%
clf; reset(gca);
surf(MUstep,Timestep,FE_Atime./FE_Dtime)
ylabel('Step Size')
xlabel('Nudging Coefficient')
zlabel('Analytic TTC')
set(gca,'YScale','log')