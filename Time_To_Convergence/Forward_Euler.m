function [M] = Forward_Euler(Matrix, Operator, N, dt, Data, MU)
M = [Matrix(:,1),zeros(size(Matrix,1),size(Matrix,2)-1)];
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

else %If Operator is a function, run this:

    if or(isequal(Data,[]),isequal(MU,[]));

        for col = 2:N-1   
        M(:,col) = M(:,col-1)+(Operator(M,col))*dt;
        end
            return;
    else

        for col = 2:N-1 
        M(:,col) = M(:,col-1)+(Operator(M,col)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
        end
            return;
    end
end

    
end