function [M] = timevo(Oper,data,x0,MU,steps)


dt = 0.0001;
Alen = size(Oper,1);


if isequal(data,[])
data = [x0;zeros(steps-1,Alen)];
for i = 2:steps
data(i,:) = abash(data,[],Oper,[],i,dt);
end
    M = data;
    return;
end

M = [x0;zeros(steps-1,Alen)];
for i = 2:steps
M(i,:) = abash(M, data, Oper, MU, i, dt);
end

end



