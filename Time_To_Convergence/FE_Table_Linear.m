% This script is just like Scheme_Time_Comparison, except it's made to loop
% through the different values to see how it all depends on different
% values. The idea is to make a table consisting of all these values to
% display in the work.
% We'll run through maybe a hundred itterations with each, which means
% running the Scheme_Time_Comparison script around 10,000 times.
% Just to test, we'll loop through the analytic and discrete first, before
% asking MATLAB to compute the time.

% NOTE: It took ~3.5 minutes to run ONLY the actual Forward Euler over 10
% different timesteps and 10 different mu values. That means if you wanted
% to plot a 100x100 graph, it would take around six hours. UHHHHH...

%% Starting the for loop
clear;clc;
% Create two 100x100 matrices for Analytic and Discrete times
% The row will represent the timestep and the column the one value of mu.
% If this is fast, we can run through changing the other values in MU.
% For the most part, initial conditions will stay the same.
EE = 10;
HH = 10;
FE_Atime = zeros(EE,HH);
FE_Dtime = zeros(EE,HH);
FE_Rtime = zeros(EE,HH);
Timestep = zeros(1,EE);
MUstep = zeros(1,HH);

x0t = [.1;.1;-.1;.1];
x0 = [-.2;-.2;-.2;.1];
                                            %x0t = rand(4,1); this was an
                                            %idea to see if it was still
                                            %mu-independent. It is.
                                            %x0 = rand(4,1);
            sx0 = x0-x0t;
            eps = 1e-16;

Operator = [-1 0.5 0.5 0.2;
            1 -2 1 1; 
            0.5 0.5 -3 -1; 
            0.5 0.5 1 -4];
Opsize = size(Operator,1);

for i = 1:EE
    for j = 0:(HH-1)
            
            dt = (1e-5)*i;
            Timestep(i) = dt;
            MU = diag([.1*j,0,0,0]);
            MUstep(j+1) = MU(1,1);
            
            
            %% Computing the Analytic Time to Convergence
            
            [P_c,D_c] = eig(Operator-MU); % Calculate the eigenvalues
            sy0_FE = (P_c^-1)*sx0; % Transform the difference into eigenvector-land
            Analytic_Time = log(eps/(norm(P_c)*norm(sy0_FE)))/log(norm(diag(exp(diag(D_c))))); % Calculate upper bound for analytic time to convergence
            
            FE_Atime(i,j+1) = Analytic_Time;
            
            %% (Forward Euler) Computing the Discrete Time to Convergence
            FE_Discrete_Time = dt/log(norm(eye(Opsize)+D_c*dt))*log(eps/((norm(P_c)*norm((P_c^-1)*sx0)))); % Calculate upper bound for discrete time to convergence
            FE_Dtime(i,j+1) = FE_Discrete_Time;
      







            %% (Forward Euler) Finding the Time to Convergence According to MATLAB

            %Analyzing the time directly:
            FE_Max_Time = max([Analytic_Time,FE_Discrete_Time]);
            
            Data1 = x0t; % Set up initial conditions for numerical scheme
            Pred1 = x0;
            
            Data1 = Forward_Euler(Data1,Operator,FE_Max_Time,dt,[],[]); % Run scheme w/o DA
            Pred1 = Forward_Euler(Pred1,Operator,FE_Max_Time,dt,Data1,MU); % Run DA using generated data
            
            Errm1 = zeros(1,size(Data1,2)); % Matrices to hold the difference and data norms 
            Data1Norm = zeros(1,size(Data1,2));
            
            for column = 1:size(Errm1,2) % Calculate norms for each step
                Errm1(column) = norm(Data1(:,column)-Pred1(:,column));
                Data1Norm(column) = norm(Data1(:,column));
            end
            
            for Num1 = 1:size(Errm1,2); % Find when MATLAB hits epsilon,
                if Errm1(Num1)<eps;
                    break;
                end
            end
            Num1 = Num1+1; % Add one because
            

            FE_Rtime(i,j+1) = Num1*dt;
               
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
surf(Timestep,MUstep,transpose(FE_Dtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Discrete TTC')

