clear;clc;clf; yalmip clear;


syst = [-1 1.5 0.5 0.2;
        1 -2 1 1; 
        0.5 0.5 -3 -1; 
        0.5 0.5 1 -4];
T = 300;
dt = 1e-5;



x0t = [5,1,-3,3];
x0m = [.2,.2,-.2,.1];


mu = [.36 0 0 0];

systda = syst;

Mt = eveuler(syst,[],x0t,[],floor(T/dt),dt); 


Mm = eveuler(syst,Mt,x0m,mu,floor(T/dt),dt);


Diff = Mm-Mt;
AbsDiff = zeros(size(Diff,1),1);
for i = 1:size(AbsDiff,1)
    AbsDiff(i) = sqrt(sum(Diff(i,:).^2));
end

semilogy(1:(size(AbsDiff,1)),AbsDiff);