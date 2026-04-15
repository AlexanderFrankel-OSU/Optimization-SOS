function [X] = abash(M, Data, Oper, MU, row, dt)
transOper = transpose(Oper);
if isequal(Data,[])
    if row <= 2
        X = M(1,:)+M(1,:)*transOper*dt;
    else
        X = M(row-1,:)+(1.5*M(row-1,:)*transOper-0.5*M(row-2,:)*transOper)*dt;
    end
return;
end
    
if row <= 2
    X = M(1,:)+(M(1,:)*transOper-MU.*(M(1,:)-Data(1,:)))*dt;
else
    X = M(row-1,:)+(1.5*M(row-1,:)*transOper-0.5*M(row-2,:)*transOper-MU.*(M(row-1,:)-Data(row-1,:)))*dt;
end


end
