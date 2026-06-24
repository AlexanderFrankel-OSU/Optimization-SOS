%% Comparison of Time to Convergence for Forward Euler and Adams-Bashforth
%This script is meant to run a comparison between the analytically
%predicted time, the Forward Euler itteration, and the actual computed
%time.
clear;clc;

%% Setup

x0t = [.1;.1;-.1;.1];
x0 = [.2;.2;-.2;.1];
sx0 = x0-x0t;
eps = 1e-16;
dt = 1e-5;

Operator = [-1 0.5 0.5 0.2;
            1 -2 1 1; 
            0.5 0.5 -3 -1; 
            0.5 0.5 1 -4];
MU = diag([.1,0,0,0]);

%% Computing the Analytic Time to Convergence

[P_c,D_c] = eig(Operator-MU); % Calculate the eigenvalues
sy0_FE = (P_c^-1)*sx0; % Transform the difference into eigenvector-land
Analytic_Time = log(eps/norm(sy0_FE))/max(real(diag(D_c))); % Calculate upper bound for analytic time to convergence

display('Analytic Time to Convergence:')
display(Analytic_Time)


%% (Forward Euler) Computing the Discrete Time to Convergence

EulerDiscOp = eye(size(Operator,1))+(Operator-MU)*dt; % Define the operator according to the numerical scheme (Forward Euler)
[P_d,D_d] = eig(EulerDiscOp); % Eigenvalues
D_v = diag(D_d); % Put eigenvalues in an array
M = max(abs(D_v)); % Take maximum magnitude eigenvalue for calculation
FE_Discrete_Time = abs(dt/log(M)*log(eps/norm(sy0_FE))); % Calculate upper bound for discrete time to convergence


display('Forward Euler Discrete Time to Convergence:')
display(FE_Discrete_Time)



%% (Forward Euler) Finding the Time to Convergence According to MATLAB

%Analyzing the time directly:
FE_Max_Time = max([Analytic_Time,FE_Discrete_Time,300]);

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

display('MATLAB says TTC for Forward Euler is:')
display(Num1*dt) % Just translating steps to real 'time'


%% (Adams-Bashforth) Computing the Discrete Time to Convergence

x1t = x0t + Operator*x0t*dt;
x1 = x0 + Operator*x0*dt; % Define the second steps using FE
% Remember sx0;
sx1 = sx0 + Operator*sx0*dt;

Opsize = size(Operator,1); % Might use everywhere

G = eye(Opsize)+(1.5*Operator-MU)*dt; % Define matrices in difference equation for AB
H = -0.5*Operator*dt;

B = [G,H;eye(Opsize),zeros(Opsize,Opsize)]; % Create block matrix
[P_ab,D_ab] = eig(B); % Diagonalize



D_abv = diag(D_ab); % Put eigenvalues in matrix
Max_EigB = max(abs(D_abv)); % Take largest eigenvalue's absolute value
z1_AB = (P_ab^-1)*[sx1;sx0]; % Create column vector by joining sx1 and sx0

AB_Discrete_Time = log((sqrt(2)*eps)/norm(z1_AB))*dt/log(Max_EigB); % Calculate discrete time to convergence for AB

display('Discrete TTC for Adams-Bashforth:')
display(AB_Discrete_Time)


%% (Adams-Bashforth) Finding the Time to Convergence According to MATLAB

Data2 = [x0t,x1t]; % Set up initial conditions for numerical scheme (requires one extra step of Euler)
Pred2 = [x0,x1];
AB_Max_Time = 10*max(Analytic_Time,AB_Discrete_Time);

Data2 = Adams_Bashforth(Data2,Operator,AB_Max_Time,dt,[],[]); % Run scheme w/o DA
Pred2 = Adams_Bashforth(Pred2,Operator,AB_Max_Time,dt,Data2,MU); % Run DA using generated data

Errm2 = zeros(1,size(Data2,2)); % Matrices to hold difference and data norms
Data2Norm = zeros(1,size(Data2,2));

for column = 1:size(Errm2,2)
    Errm2(column) = norm(Data2(:,column)-Pred2(:,column)); % Calculate norms for each step
    Data2Norm(column) = norm(Data2(:,column));
end


for Num2 = 1:size(Errm1,2); % Find when MATLAB hits epsilon
    if Errm2(Num2)<eps;
        break;
    end
end
Num2 = Num2+1; % You get it, right?

display('MATLAB says TTC for Adams-Bashforth is:')
display(Num2*dt)

semilogy((0:(size(Data2Norm,2)-1))*dt,Data2Norm)
% Maybe put a pretty table down here? Would be nice!