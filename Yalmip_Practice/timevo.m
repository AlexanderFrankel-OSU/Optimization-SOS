function [M] = timevo(Oper,data,x0,MU,steps)


dt = 0.00001;
Alen = size(Oper,1);


if data == 0
data = [x0;zeros(steps-1,Alen)];
for i = 2:steps
data(i,:) = abash(data,0,Oper,0,i,dt);
end

end

M = [x0;zeros(steps-1,Alen)];
for i = 2:steps
M(i,:) = abash(M, data, Oper, MU, i, dt);
end

end



