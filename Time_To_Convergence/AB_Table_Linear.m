
%% Starting the for loop
clear;clc;
% Create two 100x100 matrices for Analytic and Discrete times
% The row will represent the timestep and the column the one value of mu.
% If this is fast, we can run through changing the other values in MU.
% For the most part, initial conditions will stay the same.
E = 6;
HH = 40;
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
            epsilon = 1e-16;

Operator = [-1 0.5 0.5 0.2;
            1 -2 1 1; 
            0.5 0.5 -3 -1; 
            0.5 0.5 1 -4];
Opsize = size(Operator,1);

for i = 1:EE
    for j = 0:(HH-1)
            
            dt = (1e-5)*i;
            Timestep(i) = dt;
            MU = diag([.05*j,0,0,0]);
            MUstep(j+1) = MU(1,1);


            %% Computing the Analytic Time to Convergence
            
            [P_c,D_c] = eig(Operator-MU); % Calculate the eigenvalues
            sy0_FE = (P_c^-1)*sx0; % Transform the difference into eigenvector-land
            Analytic_Time = log(epsilon/(norm(P_c)*norm(sy0_FE)))/log(norm(diag(exp(diag(D_c))))); % Calculate upper bound for analytic time to convergence
            
            FE_Atime(i,j+1) = Analytic_Time;
            
           %% (Adams-Bashforth) Computing the Discrete Time to Convergence

            x1t = x0t + Operator*x0t*dt;
            x1 = x0 + Operator*x0*dt; % Define the second steps using FE
            % Remember sx0;
            sx1 = sx0 + Operator*sx0*dt;
            
            Opsize = size(Operator,1); % Might use everywhere
            
            G = eye(Opsize)+(1.5*Operator-MU)*dt; % Define matrices in difference equation for AB
            H = -0.5*Operator*dt;
            
            B = [G,H;eye(Opsize),zeros(Opsize,Opsize)]; % Create block matrix
            [P_AB,D_AB] = eig(B); % Diagonalize
            P_ABin = P_AB^(-1);
            
            
            D_abv = diag(D_AB); % Put eigenvalues in matrix
            Max_EigB = max(abs(D_abv)); % Take largest eigenvalue's absolute value
            z1_AB = (P_AB^-1)*[sx1;sx0]; % Create column vector by joining sx1 and sx0
            
            AB_Discrete_Time = dt*log(epsilon/(norm([eye(Opsize) zeros(Opsize,Opsize)]*P_ab)*norm(z1_AB)))/log(norm(D_ab)); % Calculate discrete time to convergence for AB
            
            
            AB_Dtime(i,j+1) = AB_Discrete_Time;


            %% (Adams-Bashforth) Finding the Time to Convergence According to MATLAB

            Data2 = [x0t,x1t]; % Set up initial conditions for numerical scheme (requires one extra step of Euler)
            Pred2 = [x0,x1];
            AB_Max_Time = max([Analytic_Time,AB_Discrete_Time,150]);
            
            Data2 = Adams_Bashforth(Data2,Operator,AB_Max_Time,dt,[],[]); % Run scheme w/o DA
            Pred2 = Adams_Bashforth(Pred2,Operator,AB_Max_Time,dt,Data2,MU); % Run DA using generated data
            
            Errm2 = zeros(1,size(Data2,2)); % Matrices to hold difference and data norms
            Data2Norm = zeros(1,size(Data2,2));
            
            for column = 1:size(Errm2,2)
                Errm2(column) = norm(Data2(:,column)-Pred2(:,column)); % Calculate norms for each step
                Data2Norm(column) = norm(Data2(:,column));
            end
            
            
            for Num2 = 1:size(Errm2,2); % Find when MATLAB hits epsilon
                if Errm2(Num2)<epsilon;
                    break;
                end
            end
            Num2 = Num2+1; % You get it, right?

            AB_Rtime(i,j+1) = Num2*dt;
    end
end
%% Display
%display(1+log(epsilon/norm(Lcpp))/log(Max_EigB))
%display(1+log(epsilon/norm(Lc13))/log(Max_EigB))
%{
AnDidiff = abs(AB_Atime-AB_Dtime);
ABlogAnDidiff = log(AnDidiff);

surf(Timestep,MUstep,transpose(AB_Dtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Discrete TTC')

surf(Timestep,MUstep,transpose(AB_Atime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Analytic TTC')

surf(Timestep,MUstep,transpose(ABlogAnDidiff))
xlabel('Timestep')
ylabel('MU value')
zlabel('Log of difference in TTC')

surf(Timestep,MUstep,transpose(AB_Rtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Real TTC')
%}

surf(Timestep,MUstep,transpose(AB_Dtime))
xlabel('Timestep')
ylabel('MU Value')
zlabel('Discrete TTC')
