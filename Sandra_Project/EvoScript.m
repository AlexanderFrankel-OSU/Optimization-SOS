function EvoScript(Ux,Uy,Uz)

% If there's a better way to reassign the variables ux,uy,uz without
% creating new variables, I would be very excited to know. I found that
% this works, however.

load("Sandra_Vars.mat");
ux=Ux;uy=Uy;uz=Uz;
clear("Ux","Uy","Uz");
save("Sandra_Vars.mat");
load("Sandra_Vars.mat");

daPos=initmatrix(daPos);
Pos=initmatrix(Pos);

Pos = adbash(Pos);
daPos = DAadbash(daPos,Pos);

% Create a difference variable  and use it to define the length of the
% difference vector
diffP = daPos - Pos; 
absdiff = (diffP(:,1).^2+diffP(:,2).^2+diffP(:,3).^2).^0.5;

% Plot the absolute difference evolution over the number of itterations.
semilogy(absdiff);

daPos = initmatrix(daPos);
Pos = initmatrix(Pos);

end