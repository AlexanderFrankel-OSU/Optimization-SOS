function [M] = Adams_Bashforth(Matrix, Operator, N, dt, Data, MU)

M = [Matrix(:,1:2),zeros(size(Matrix,1),size(Matrix,2)-2)];

if isequal(M(:,2),zeros(size(M,1),1));
    display('Enter initial state')
   return;
end

%If Operator is a matrix, run this:

if isequal(class(Operator),'double');

    if or(isequal(Data,[]),isequal(MU,[]));
        
        for col = 3:N+1   
        M(:,col) = M(:,col-1)+1.5*Operator*(M(:,col-1)-0.5*M(:,col-2))*dt;
        end
            return;
    else

        for col = 3:N+1 
        M(:,col) = M(:,col-1)+(Operator*(1.5*M(:,col-1)-0.5*M(:,col-2))-MU*(M(:,col-1)-Data(:,col-1)))*dt;
        end
        return;
    end

else %If Operator is a function, run this:

    if or(isequal(Data,[]),isequal(MU,[]));

        for col = 3:N+1   
        M(:,col) = M(:,col-1)+(1.5*Operator(M,col-1)-0.5*Operator(M,col-2))*dt;
        end
            return;
    else

        for col = 3:N+1 
        M(:,col) = M(:,col-1)+(1.5*Operator(M,col-1)-0.5*Operator(M,col-2)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
        end
            return;
    end
end

    
end


