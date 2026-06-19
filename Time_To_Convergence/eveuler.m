function [M] = eveuler(Oper,Data,x0,MU,steps,dt)

Oplen = size(Oper,1); 

if or(isequal(Data,[]),isequal(MU,[]))

Data = [x0;zeros(steps-1,Oplen)]; 
    for i = 2:steps
    Data(i,:) = feuler(Data,[],Oper,[],i,dt);
    end
        M = Data; 
        return;
end


M = [x0;zeros(steps-1,Oplen)];

for i = 2:steps
    M(i,:) = feuler(M, Data, Oper, MU, i, dt);
end

end



