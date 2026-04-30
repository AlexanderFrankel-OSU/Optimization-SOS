clear;clc;clf;
yalmip clear

% Define the update operator as a matrix so that abash can evaluate the
% timestep using dx/dt = Ax
syst = [-1 1.5 0.5 0.2;
        1 -2 1 1; 
        0.5 0.5 -3 -1; 
        0.5 0.5 1 -4];
T = 200;
dt = 1e-5;


% In this script the 't' stands for true, and the 'm' stands for modeled,
% so the initial data here shows we're starting somewhere fairly different
% from where we should be:
x0t = [5,1,-3,3];
x0m = [.2,.2,-.2,.1];

% As the IntroAlgorithm script seemed to suggest, the second mu value
% approached zero for MOSEK, so I set only the first one close to its
% minimum:
mu = [.36 0 0 0];

systda = syst;%-diag(mu)

% Note that the data assimilation system below was what I was using before; 
% I just decided to switch to the diagonal because it was cleaner, and can 
% probably be applied in the future.

% [-(1+mu(1)) 1.5 0.5 0.2; 
%        1 -(2+mu(2)) 1 1; 
%        0.5 0.5 -3 -1; 
%        0.5 0.5 1 -4];

% Here's where I make use of my (likely faulty) time evolution function,
% which uses the abash function.
% Entering [] where the data goes (second entry) or [] where the mu vector
% goes (fourth entry) will cause the time evolution function to generate
% the data without data assimilation.
% The operator matrix goes in the first entry, initial values in the third,
% and number of abash itterations in the last.

Mt = timevo(syst,[],x0t,[],floor(T/dt),dt); 
% I'll add shorter comments in the timevo and abash scripts.

% Once the "data" is generated from the above line, the model takes in the
% data "Mt" and uses the mu values, starting at different initial values.
Mm = timevo(syst,Mt,x0m,mu,floor(T/dt),dt);

% After that, we take the difference between data and model and calculate 
% the absolute difference in quadrature, and plot it using semilogy.
Diff = Mm-Mt;
AbsDiff = zeros(size(Diff,1),1);
for i = 1:size(AbsDiff,1)
    AbsDiff(i) = sqrt(sum(Diff(i,:).^2));
end

semilogy(1:(size(AbsDiff,1)),AbsDiff);

%Turn on if you want to see what's happening to the "true system":
%Maxim = zeros(size(Mt,1),1);
%for i = 1:size(Mt,1)
%    Maxim(i) = sqrt(sum(Mt(i,:).^2));
%end
%plot(1:size(Maxim,1),Mt)
