clear;clc;
%This script is meant to run a comparison between the analytically
%predicted time, the Forward Euler itteration, and the actual computed
%time.


x0t = [5;1;-3;3];
x0 = [.2;.2;-.2;.1];
sx0 = x0-x0t;

eps = 1e-16;

Operator = [-1 1.5 0.5 0.2;
            1 -2 1 1; 
            0.5 0.5 -3 -1; 
            0.5 0.5 1 -4];
MU = diag([.36,0,0,0]);
%Computing the Analytic Time to Convergence
[P,D] = eig(Operator-MU);
sy0 = (P^-1)*sx0;
dt = 1e-5;

Analytic_Time = log(eps/norm(sy0))/max(real(diag(D)));



%Computing the Discrete Time to Convergence
EulerDiscOp = eye(size(Operator,1))-(Operator-MU)*dt;
[P_d,D_d] = eig(EulerDiscOp);
D_v = diag(D_d);
M = max(abs(D_v));

Discrete_Time = abs(dt/log(M)*log(eps/norm(sy0)));

display('Analytic time to convergence is:')
display(Analytic_Time)
%display('The norm of sx is:')
%display(norm(P*exp(D*Analytic_Time)*(P^-1)*sx0))
%display('The norm of sy is:')
%display(norm(exp(D*Analytic_Time)*sy0))
display('Discrete time to convergence is:')
display(Discrete_Time)

%N_an = ceil(Analytic_Time/dt);
%N_dsc = ceil(Discrete_Time/dt);

%Matrix_tan = [x0t,zeros(size(x0t,1),N_an)];
%Matrix_an = [x0,zeros(size(x0,1),N_an)];
%Data_an = Forward_Euler(Matrix_tan, Operator, N_an, dt, [], []);
%Results = Forward_Euler(Matrix_an, Operator, N_an, dt, Data, MU);


%Matrix_tdsc = [x0t,zeros(size(x0t,1),N_dsc)];
%Matrix_dsc = [x0,zeros(size(x0t,1),N_dsc)];
%Data_dsc = Forward_Euler(Matrix_tdsc, Operator, N_dsc, dt, [], []);
%Results_dsc = Forward_Euler(Matrix_dsc, Operator, N_dsc, dt, Data_dsc, MU);

%Analyzing the time directly:
Max_Time = max(Analytic_Time,Discrete_Time);
Max_N = ceil(Max_Time/dt);

Data1 = [x0t,zeros(size(Operator,1),Max_N)];
Pred1 = [x0,zeros(size(Operator,1),Max_N)];


Data1 = Forward_Euler(Data1,Operator,Max_N,dt,[],[]);
Pred1 = Forward_Euler(Pred1,Operator,Max_N,dt,Data1,MU);
Errm1 = zeros(1,size(Data1,2));
for i = 1:size(Errm1,2)
Errm1(i) = norm(Data1(:,i)-Pred1(:,i));
end


Data2 = [x0t,x0t+Operator*x0t*dt,zeros(size(Operator,1),Max_N-1)];
Pred2 = [x0,x0+Operator*x0*dt,zeros(size(Operator,1),Max_N)-1];
Data2 = Adams_Bashforth(Data2,Operator,Max_N,dt,[],[]);
Pred2 = Adams_Bashforth(Pred2,Operator,Max_N,dt,Data1,MU);
Errm2 = zeros(1,size(Data1,2));
for i =1:size(Errm2,2)
    Errm2(i)=norm(Data2(:,i)-Pred2(:,i));
end

