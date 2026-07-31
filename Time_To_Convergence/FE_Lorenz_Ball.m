% This script runs simulations of the Lorenz '63 System with data
% assimilation using minimal nudging coefficients.

%%
clear;clc;

sig = 10;
beta = 8/3;
rho = 28;
dt = 1e-5;
s01 = [2;-3;1.5];
s02 = [0;0;0];

xmu = 10.9795971282064; % Minimal nudging coefficients
ymu = 14.6854347255662;
zmu = 0.00127537481573532;

MUcmb = [10.97960 14.68543 0.001275376; ...
         10.97960 14.68543 0; ...
         xmu 0   zmu; ...
         0   ymu zmu; ...
         319.0741 0   0; ...
         0   31.90741 0; ...
         0   0   zmu];
CMBS = size(MUcmb,1);



dSdt = @(s) [-sig*(s(1)-s(2)); rho*s(1)-s(2)-s(1)*s(3); s(1)*s(2)-beta*s(3)];
%{
%% Running a Single Simulation
max_time = 200;

CMBO = 5;
    MU = diag(MUcmb(CMBO,:));
    LData = FrdEulCntrlTD(s01,dSdt,max_time,dt);
    LPred = FrdEulCntrlTD(s02,dSdt,max_time,dt,LData,MU);
    LDiff = LPred-LData;
    TTC = TTCfinder(LDiff,dt);
    TTC_max = max(TTC);


%% Plot?
timeline = (0:size(LDiff,2)-1)*dt;
normLDiff = vecnorm(LDiff,2,1);
semilogy(timeline,normLDiff)
% plot(timeline,LDiff)
% plot3(LDiff(1,:),LDiff(2,:),LDiff(3,:))
% plot3(LData(1,:),LData(2,:),LData(3,:))
%}



%% Figure-Making Code

NumTS = 20; % Number of timesteps
stepend = 1e-1;
rcoff = (stepend/1e-5)^(1/(NumTS-1)); % Geometrically distributed between 1e-5 and 0.1.
Timestep = transpose((1e-5)*rcoff.^(0:NumTS-1)); % Range of timesteps to compute over
EE = length(Timestep);



HH = 20; % However many mu values for x
VV = 1; % Number of mu values for y
MUstep = cell(HH,VV);

% If running only x, or only y, set the corresponding value to 1. 
% Don't forget to change the value of MUstep.

% % for both x and y
% MUline = zeros(1,HH);
% for JJ = 1:HH
%     for KK = 1:VV
%         MUstep{JJ,KK} = [MUcmb(2,1)/HH*JJ MUcmb(2,2)/VV*KK 0];
%     end
% end

% For x only
MUline = zeros(1,HH);
for JJ = 1:HH
        MUstep{JJ} = [20 0 0]/HH*JJ;%MUcmb(5,:)/HH*JJ;
        MUline(JJ) = MUstep{JJ}(1);
end

% % For y only
% MUline = zeros(1,VV);
% for KK = 1:VV
%         MUstep{KK} = MUcmb(6,:)/VV*KK;
%         MUline(KK) = MUstep{KK}(2);
% end

%%
TIME = 400;
[FE_Rtime Convgs CFL] = deal(cell(1,EE)); % Set up the multiindex cell arrays!
for II = 1:EE
    [FE_Rtime{II} Convgs{II}] = deal(zeros(HH,VV));
    CFL{II} = ones(HH,VV);

end

%%

    for II = 1:EE
        for JJ = 1:HH
            for KK = 1:VV
                if all(MUstep{JJ}<2/Timestep(II)) % Check CFL condition. Don't run unless it isn't violated.
                    dt = Timestep(II);
                    MU = diag(MUstep{JJ,KK});
                    part1 = FrdEulCntrlTD(s01,dSdt,TIME,dt);
                    part2 = FrdEulCntrlTD(s02,dSdt,TIME,dt,part1,MU);
                    pdiff = part2-part1;
                    FE_Rtime{II}(JJ,KK) = max(TTCfinder(pdiff,dt));
                    % Shouldn't converge at the end, or blow up.
                    if ~(FE_Rtime{II}(JJ,KK)<TIME) || ~(all(abs(part1(:,end))<1e+10)) || ~(all(abs(part2(:,end))<1e+10))
                        FE_Rtime{II}(JJ,KK) = NaN;
                    else
                        Convgs{II}(JJ,KK) = 1;
                    end

                else
                    FE_Rtime{II}(JJ,KK) = NaN;
                    CFL{II}(JJ,KK) = 0;
                end
            end
        end
    end


display('We made it!!!')


%% Plot!!
% Only run the next for loop if doing only x or only y.
% Convert the cells into matrices, like they were always meant to be.
%%
[FE_RtimeMAT ConvgsMAT CFLMAT] = deal(zeros(EE,max(HH,VV)));
for II = 1:EE
    for JJ = 1:HH
        for KK = 1:VV
            FE_RtimeMAT(II,max(JJ,KK)) = FE_Rtime{II}(JJ,KK);
            ConvgsMAT(II,max(JJ,KK)) = Convgs{II}(JJ,KK);
            CFLMAT(II,max(JJ,KK)) = CFL{II}(JJ,KK);
        end
    end
end


clf;
surf(MUline,Timestep,FE_RtimeMAT) % Clears the outliers.
ylabel('Step Size')
xlabel('Nudging Coefficient')
zlabel('Numerical TTC')
set(gca,'YScale','log')


%% Not to be run. Used originally to calculate minimal MU values.
%{
% % This script is meant to calculate the minimal nudging coefficients for
% the Lorenz system when we know that the system will at some point be
% within an absorbing ball (because the Lorenz '63 system is dissipative).
clear; clc; yalmip clear;
negdefconst = 1e-9;

% Define the decision variables
mu = sdpvar(3,1);
X = sdpvar(3,1);
Xtrue = sdpvar(3,1);
% [-2.3365;3.263;1.7192];
optvar = mu;


% Define related constants
sig = 10;
beta = 8/3;
rho = 28;

% Here's all the YALMIP stuff

sym_matrix = [-(sig+mu(1)) (sig+rho-Xtrue(3))/2 Xtrue(2)/2;
              (sig+rho-Xtrue(3))/2 -(1+mu(2)) 0;
              Xtrue(2)/2 0 -(beta+mu(3))];
dUdt = X'*sym_matrix*X

neg_dUdt_sos = -dUdt-negdefconst*dot(X,X); % Negative definiteness constraint

obj_func = dot(mu,mu);
% Introduce the constraints
% Hmm, not working so well: , X(1)^2 <= 2.3365^2, X(2)^2 <= 3.263^2, X(3) >= 0, X(3) <= 1.7192
constr = [sos(neg_dUdt_sos), mu(:) >= 0, Xtrue(1)^2 <= 2.3365^2, Xtrue(2)^2 <= 3.263^2, Xtrue(3) >= 0, Xtrue(3) <= 1.7192];

options = sdpsettings('solver', 'mosek');

[sol,v,M,res] = solvesos(constr,obj_func,options,optvar);
format long g
value(mu)
%}