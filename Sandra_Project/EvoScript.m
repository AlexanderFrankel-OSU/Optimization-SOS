load("Sandra_Vars.mat")

daPos=initmatrix(daPos);
Pos=initmatrix(Pos);

Pos = adbash(Pos);
daPos = DAadbash(daPos,Pos);

diffP = daPos - Pos; 
absdiff = (diffP(:,1).^2+diffP(:,2).^2+diffP(:,3).^2).^0.5;

semilogy(absdiff);

daPos = initmatrix(daPos);
Pos = initmatrix(Pos);

