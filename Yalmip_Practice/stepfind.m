clear;clc;clf;yalmip clear;
% Basic idea:  dxdt = x*A --> dxdt = xPDP^(-1) --> (dxdt*P)=xPD --> dydt = y*D
% d(Sy)dt = Sy*(D-diag(mu))
syst = [-1 1.5 0.5 0.2;
        1 -2 1 1; 
        0.5 0.5 -3 -1; 
        0.5 0.5 1 -4];

mu = [.36 0 0 0];
diagmu = diag(mu);
x0t = [5,1,-3,3];
x0m = [.2,.2,-.2,.1];
Oper = syst;

[Evecs,Evals] = eig(Oper);
SOper = Evals-diagmu;

y01 = x0t*(transpose(Evecs))^(-1);
y02 = x0m*(transpose(Evecs))^(-1);
Sy0 = y02-y01;

T=100;
dt = 1e-5;

% The operator is a diagonal matrix... It does not change under transpose!
SyEvo = timevo(SOper,[],Sy0,[],floor(T/dt),dt);
SxEvo = SyEvo*transpose(Evecs);
Diff0 = zeros(size(SxEvo,1),1);
for i = 1:size(SxEvo,1)
    Diff0(i) = sqrt(sum(real(SxEvo(i,:)).^2));
end
LogDiff0 = log10(Diff0);
Timeplogerr = [(1:size(LogDiff0,1))'*dt,LogDiff0]; 

trunc = Timeplogerr(floor(.1*T/dt):end,:);

mxpb = polyfit(trunc(:,1),trunc(:,2),1);
suffT = (-12-mxpb(2))/mxpb(1);
plot((1:size(LogDiff0,1))*dt,LogDiff0)
% plot((1:size(LogDiff0,1))*dt,mxpb(1)*((1:size(LogDiff0,1))*dt)+mxpb(2))

display('Based on the system, nudging coefficients, and timestep, a reasonable expectation for the amount of time needed is:')
display(suffT)

% d(Sy)dt = Sy*(D-diag(mu));





