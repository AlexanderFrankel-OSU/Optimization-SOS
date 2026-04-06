function [X] = abash(M, Data, Oper, MU, row, dt)
if Data == 0
    if row <= 2
        X = (M(1,:)'+(Oper*(M(1,:)'))*dt)';
    else
        X = (M(row-1,:)'+(1.5*Oper*(M(row-1,:)')-0.5*Oper*(M(row-2,:)'))*dt)';
    end
else
    if row <= 2
        X = (M(1,:)'+(Oper*(M(1,:)')-MU.*(M(1,:)'-Data(1,:)'))*dt)';
    else
        X = (M(row-1,:)'+(1.5*Oper*(M(row-1,:)')-0.5*Oper*(M(row-2,:)')-MU.*(M(row-1,:)'-Data(row-1,:)'))*dt)';
    end
end

end
