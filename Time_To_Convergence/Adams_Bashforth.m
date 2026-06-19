function [X] = Adams_Bashforth(M, Data, Oper, MU, row, dt)
% Runs one Adams-Bashforth step at a specified row.

% Because of how I decided to store data in a matrix, we need to transpose
% the operator which is written to generate data in column vector format:
transOper = transpose(Oper);
% My idea was to take d/dt([x1,x2,...,xn]^T)=Oper*[x1,x2,...,xn]^T 
% and transpose the equation so that we have:
% d/dt([x1,x2,...,xn]) = [x1,x2,...,xn]*(Oper^T), and can just multiply on
% the left by the row vector.

% Here, check if the Data or MU entries are excluded and run without data
% assimilation if so.
if or(isequal(Data,[]),isequal(MU,[]))
    % Check if the row is 2 or less because Adams-Bashforth requires two
    % rows before it, so if not, run one step of Euler.
    if row <= 2
        X = M(1,:)+M(1,:)*transOper*dt;
    else
    % Otherwise, run Adams-Bashforth
        X = M(row-1,:)+(1.5*M(row-1,:)*transOper-0.5*M(row-2,:)*transOper)*dt;
    end
return;
end

% If the Data and MU are full, then run with data assimilation:
if row <= 2
    X = M(1,:)+(M(1,:)*transOper-MU.*(M(1,:)-Data(1,:)))*dt;
else
    X = M(row-1,:)+(1.5*M(row-1,:)*transOper-0.5*M(row-2,:)*transOper-MU.*(M(row-1,:)-Data(row-1,:)))*dt;
end


end
