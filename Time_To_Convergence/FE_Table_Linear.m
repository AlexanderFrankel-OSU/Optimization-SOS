% This script loops through the different nudging coefficients and timesteps
% to see how it all depends on different values. The idea is to make a table 
% consisting of all these values to display in the work.

% NOTE: It took ~3.5 minutes to run ONLY the actual Forward Euler over 10
% different timesteps and 10 different mu values. If you wanted
% to plot a 100x100 graph, it would take around six hours.

%% Setting up the necessary initial conditions


clear;clc;


x0t = [.1;.1;-.1;.1];
x0 = [-.2;-.2;-.2;.1];

sx0 = x0-x0t;
dt = 1e-5;
epsilon = 1e-16;
% The algorithm gives the value 0.823543613122324 for the minimal nudging
% coefficient. However, in practice, the minimum should be exactly 0.62.
% We'll run starting at 0.63
mumin = 0.823543613122324;
trumin = 0.62;

Operator = [-1 2 0.5 0.2;
             1 -2 1 1;
             0.5 0.5 -3 -1;
             0.5 0.5 1 -4];


MU = diag([trumin+.05 0 0 0]);

display('Note the matrix here causes uncontrolled growth in variables.')
display('The eigenvalues are:')
display(eig(Operator))
display('The minimal nudging coefficient is:')
display(mumin(1))
Opsize = size(Operator,1);

% The code below can be ignored if you just want to generate the table.
% If that's your goal, skip to the section "Begin long for loop"
%{

%% Computing the Analytic Time to Convergence

[P_c,D_c] = eig(Operator-MU); % Calculate the eigenvalues and eigenvectors
Analytic_Time = log(epsilon/(norm(P_c)*norm((P_c^-1)*sx0)))/log(norm(expm(D_c))); % Calculate upper bound for analytic time to convergence
% Note, the two extra diag()'s in the log is because MATLAB calculates exp()
% element-wise, which creates ones when you take exp(Diagonal_Matrix),
% which interferes with the operator norm.
display('Analytic Time to Convergence:')
display(Analytic_Time)


%% (Forward Euler) Computing the Discrete Time to Convergence

FE_Discrete_Time = dt/log(norm(eye(Opsize)+D_c*dt))*log(epsilon/((norm(P_c)*norm((P_c^-1)*sx0)))); % Calculate upper bound for discrete time to convergence


display('Forward Euler Discrete Time to Convergence:')
display(FE_Discrete_Time)

%% Run one round of Forward Euler
FE_Max_Time = max([Analytic_Time,FE_Discrete_Time]); % Arbitrarily runs for twice the predicted bound, just to catch any edge cases

Data1 = Forward_Euler(x0t,Operator,FE_Max_Time,dt,[],[]); % Run scheme w/o DA
Pred1 = Forward_Euler(x0,Operator,FE_Max_Time,dt,Data1,MU); % Run DA using generated data

%% Better TTC Catcher

lindiff = Pred1-Data1;
TTC = zeros(Opsize,1);
if anynan(lindiff(:,end))
     notNAN = zeros(Opsize,1);
     startNAN = size(lindiff,2);
     for row = 1:Opsize
        notNAN(row,1) = find(~isnan(lindiff(row,1:startNAN)),1,"last");
        startNAN = notNAN(row,1);
     end
     notNANmin = min(notNAN);
     for row = 1:Opsize
        TTC(row,1) = (find(lindiff(row,1:notNANmin),1,"last")+1)*dt;
     end
     TTC_max = max(TTC);
else
     for row = 1:Opsize
         TTC(row,1) = (find(lindiff(row,:),1,"last")+1)*dt;
     end
     TTC_max = max(TTC);
end

timeline = (0:size(lindiff,2)-1)*dt;
%}

%% Begin long for loop

% Create two matrices for Analytic and Discrete times
% The row will represent the timestep and the column the one value of mu.
% If this is fast, we can run through changing the other values in MU.
% For the most part, initial conditions will stay the same.

rcoff = (1e+5)^(1/10);
Timestep = (1e-5)*rcoff.^(0:10);
EE = size(Timestep,2);
HH = 50;
barr = 0.02;
% MUstep = trumin+(barr:(mumin-trumin-barr)/(HH-1):mumin-trumin); % Only to plot between the true minimum and the MOSEK minimum.
mumax = 10; % Otherwise run this for the larger, less-interesting plot.
MUstep = mumin+(0:mumax/(HH-1):mumax);

FE_Atime = zeros(EE,HH);
FE_Dtime = zeros(EE,HH);
FE_Rtime = zeros(EE,HH);
TimetoRunFE = zeros(EE,HH);



for II = 1:EE
    for JJ = 1:HH % Start at zero so that we see no convergence occurs if mu is zero
            
            dt = Timestep(1,II);
            MU = diag([MUstep(1,JJ),0,0,0]);
            
            
            %% Computing the Analytic Time to Convergence
            
            [P_c,D_c] = eig(Operator-MU); % Calculate the eigenvalues
            sy0_FE = (P_c^-1)*sx0; % Transform the difference/ into eigenvector-land
            Analytic_Time = log(epsilon/(norm(P_c)*norm(sy0_FE)))/log(norm(expm(D_c))); % Calculate upper bound for analytic time to convergence
            if Analytic_Time < 0
                FE_Atime(II,JJ) = NaN;
            else
                FE_Atime(II,JJ) = Analytic_Time;
            end
            %% (Forward Euler) Computing the Discrete Time to Convergence
            FE_Discrete_Time = dt/log(norm(eye(Opsize)+D_c*dt))*log(epsilon/((norm(P_c)*norm((P_c^-1)*sx0)))); % Calculate upper bound for discrete time to convergence
            FE_Dtime(II,JJ) = FE_Discrete_Time;
            
            if FE_Discrete_Time < 0
                FE_Dtime(II,JJ) = NaN;
            else
                FE_Dtime(II,JJ) = FE_Discrete_Time;
            end

            %% (Forward Euler) Finding the Time to Convergence According to MATLAB
            
            %Analyzing the time directly:
            FE_Max_Time = max([Analytic_Time,FE_Discrete_Time]);
            
            Data1 = x0t; % Set up initial conditions for numerical scheme
            Pred1 = x0;

            tic;
            Data1 = Forward_Euler(Data1,Operator,FE_Max_Time,dt,[],[]); % Run scheme w/o DA
            Pred1 = Forward_Euler(Pred1,Operator,FE_Max_Time,dt,Data1,MU); % Run DA using generated data
            TimetoRunFE(II,JJ) = toc;

            lindiff = Pred1-Data1;
            TTC = zeros(Opsize,1);
            if anynan(lindiff(:,end))
                 notNAN = zeros(Opsize,1);
                 startNAN = size(lindiff,2);
                 for row = 1:Opsize
                    notNAN(row,1) = find(~isnan(lindiff(row,1:startNAN)),1,"last");
                    startNAN = notNAN(row,1);
                 end
                 notNANmin = min(notNAN);
                 for row = 1:Opsize
                    TTC(row,1) = (find(lindiff(row,1:notNANmin),1,"last")+1)*dt;
                 end
            else
                for row = 1:Opsize
                    TTC(row,1) = (find(lindiff(row,:),1,"last")+1)*dt;
                end
            end

            FE_Rtime(II,JJ) = max(TTC);  
            
    end
end

%{
AnDidiff = abs(FE_Atime-FE_Dtime);
logAnDidiff = log(AnDidiff);

surf(Timestep,MUstep,transpose(FE_Dtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Discrete TTC')

surf(Timestep,MUstep,transpose(FE_Atime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Analytic TTC')

surf(Timestep,MUstep,transpose(logAnDidiff))
xlabel('Timestep')
ylabel('MU value')
zlabel('Log of difference in TTC')
%}
%%
surf(Timestep,MUstep,transpose(FE_Rtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Real TTC')

%%

% surf(Timestep,MUstep,transpose(FE_Dtime))
% xlabel('Timestep')
% ylabel('MU Value')
% zlabel('Discrete TTC')


%%
%{
surf(Timestep,MUstep,transpose(abs(FE_Rtime-FE_Dtime)))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Log of difference between discrete TTC and Actual')
%}

%%
%{
surf(Timestep,MUstep,transpose(FE_Atime./FE_Rtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Ratio of analytic time to convergence to numerical')
%}

%% Code to Produce the Minimal Nudging Coefficients (Only Run if Changed)
%{
clear;clc;
yalmip clear

% Define sdpvars x and mu, and set mu to be the optimization variable
% (we're not tweaking the p_s or q_s [which are the four entries of x]).
x = sdpvar(4,1);
mu = sdpvar(1,1);
optvar = mu;
% Create a small epsilon to ensure definiteness!
epsilon = 1e-8;

% Define Q as a symmetric matrix (by turning a triangular sum of quadratic
% terms into the form z'*Q*z from the system).
% KUE =       [-1 2 0.5 0.2;
%               1 -2 1 1;
%               0.5 0.5 -3 -1;
%               0.5 0.5 1 -4]

KUE = [-1 2 0.5 0.2;
             1 -2 1 1;
             0.5 0.5 -3 -1;
             0.5 0.5 1 -4];


Q = (KUE+KUE')/2-diag([mu(1),0,0,0]);

dUdt = x'*Q*x;
%Subtract a small fraction of the sum of squares to ensure definiteness
negdUdtsos = -dUdt;

% Define other constraints, such as that
objfunc = dot(mu,mu);

% Define the constraints on the optimization. The mu_s need to be strictly
% positive, and Q should be negative semi-definite (which is equivalent to
% the expression for dUdt being negative semi-definite, which MOSEK can
% handle).
constr = [sos(negdUdtsos)];

options = sdpsettings('solver','mosek');

[sol,v,M,res] = solvesos(constr,objfunc,options,optvar);

eig(KUE)

format long g;
value(mu)
format default;
%}