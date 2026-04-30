function [X] = feuler(M, Data, Oper, MU, row, dt)

transOper = transpose(Oper);

if or(isequal(Data,[]),isequal(MU,[]))
        X = M(row-1,:)+(M(row-1,:)*transOper)*dt;
        return;
end

    X = M(row-1,:)+(M(row-1,:)*transOper-MU.*(M(row-1,:)-Data(row-1,:)))*dt;

end
