clear;clc;
yalmip clear

syst = [-1 1.5 0.5 0.2;
        1 -2 1 1; 
        0.5 0.5 -3 -1; 
        0.5 0.5 1 -4];

x0t = [5,1,-3,3];
x0m = [.2,.2,-.2,.1];

mu = [.36 0 0 0];

systda = syst-diag(mu);
% [-(1+mu(1)) 1.5 0.5 0.2;
%        1 -(2+mu(2)) 1 1; 
%        0.5 0.5 -3 -1; 
%        0.5 0.5 1 -4];

Mt = timevo(syst,[],x0t,[],1e+6);
Mm = timevo(systda,Mt,x0m,mu,1e+6);

Diff = Mm-Mt;
Diffsqr = Diff.^2;
AbsDiff = zeros(size(Diffsqr,1),1);
for i = 1:size(AbsDiff,1)
    AbsDiff(i) = sqrt(sum(Diffsqr(i,:)));
end



semilogy(1:1e+6,AbsDiff);