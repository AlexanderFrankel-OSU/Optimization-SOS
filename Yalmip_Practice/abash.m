function [X] = abash(M, Data, Oper, MU, row, dt)

if isequal(Data,[])
    if row <= 2
        X = M(1,:)+M(1,:)*Oper*dt;
    else
        X = M(row-1,:)+(1.5*M(row-1,:)*Oper-0.5*M(row-2,:)*Oper)*dt;
    end
return;
end
    
if row <= 2
    X = M(1,:)+(M(1,:)*Oper-MU.*(M(1,:)-Data(1,:)))*dt;
else
    X = M(row-1,:)+(1.5*M(row-1,:)*Oper-0.5*M(row-2,:)*Oper-MU.*(M(row-1,:)-Data(row-1,:)))*dt;
end


end
