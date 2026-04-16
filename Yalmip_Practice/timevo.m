function [M] = timevo(Oper,Data,x0,MU,steps)
% Returns the calculated time evolution matrix.

% This time step is arbitrary, and can be made smaller; I can also change
% the scripts to allow it as an optional last parameter.
dt = 1e-5; 
% Decreasing to 1e-4 makes it converge faster, of course, but I worry about violating the CFL condition.
Oplen = size(Oper,1); 
% Taking the length of the operator tells the script how many columns to
% expect (how many variables).

% This if statement is meant to run if we disinclude the data or mus.
% It then itterates through the number of inputted steps with abash.
if or(isequal(Data,[]),isequal(MU,[]))
% Initialize the reused variable Data with the starting position
Data = [x0;zeros(steps-1,Oplen)]; 
% Then run abash on all the following rows of the Data matrix:
for i = 2:steps
Data(i,:) = abash(Data,[],Oper,[],i,dt);
end
    M = Data; % If this if statement runs, it returns the matrix M.
    return;
end

% If both slots for data and mu are full, then run abash on the output
% matrix M, and give it the initial data x0.
M = [x0;zeros(steps-1,Oplen)];

for i = 2:steps
% Take in the values provided and compute one round of Adams-Bashforth
% until we reach the final row:
M(i,:) = abash(M, Data, Oper, MU, i, dt);
end
% At this point, whatever M is should be full and is returned.
end



