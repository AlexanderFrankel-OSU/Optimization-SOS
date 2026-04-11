clear;clc;
yalmip clear

syst = [-1 1.5 0.5 0.2;
        1 -2 1 1; 
        0.5 0.5 -3 -1; 
        0.5 0.5 1 -4];

x0t = [2,2,-2,1];
x0m = [.1,-.2,-.5,.2];

mu = [10;2;0;0];

systda = [-(1+mu(1)) 1.5 0.5 0.2;
    1 -(2+mu(2)) 1 1; 
    0.5 0.5 -3 -1; 
    0.5 0.5 1 -4];

Mt = timevo(syst,[],x0t,[],8e+5);
Mm = timevo(systda,Mt,x0m,mu,8e+5);

Diff = Mm-Mt;
Diffsqr = Diff.^2;
Diffsqr = sum(Diffsqr')';


semilogy(Diffsqr);