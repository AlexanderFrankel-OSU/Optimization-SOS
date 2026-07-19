function [M] = Forward_Euler(InitCons, Operator, Time, dt, Data, MU)
N = ceil(Time/dt);
if size(InitCons,1)*N*8 > 1.41e+10
    M = NaN;
    return;
end
M = [InitCons(:,1),zeros(size(InitCons,1),N)];
if isequal(M(:,1),zeros(size(M,1),1));
    display('Enter initial state')
   return;
end
%If Operator is a matrix, run this:
if isequal(class(Operator),'double');

    if or(isequal(Data,[]),isequal(MU,[]));

        for col = 2:N+1   
            M(:,col) = M(:,col-1)+(Operator*M(:,col-1))*dt;
        end
            return;
    else

        for col = 2:N+1 
            M(:,col) = M(:,col-1)+(Operator*M(:,col-1)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
        end
            return;
    end

elseif isequal(class(Operator),'function_handle') %If Operator is a function, run this:

    if or(isequal(Data,[]),isequal(MU,[]));
        if nargin(Operator)<2
            for col = 2:N+1   
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1)))*dt;
            end
        elseif nargin(Operator)==2
            for col = 2:N+1   
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1),Data(:,col-1)))*dt;
            end
        else 
                display('Check the number of variables in your differential operator.')
            return;
        end
                return;
    else
        if nargin(Operator)<2
            for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1))-MU*(M(:,col-1)-Data(:,col-1)))*dt;
            end
        elseif nargin(Operator)==2
            for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1),Data(:,col-1))-MU*(M(:,col-1)-Data(:,col-1)))*dt;
            end
        else
                display('Check the number of variables in your differential operator.')
        end
            return;
    end
end

    
end