for i=3:500000
Pos(i,1)=Pos(i-1,1)+1.5*sig*(Pos(i-1,2)-Pos(i-1,1))*dt-0.5*sig*(Pos(i-2,2)-Pos(i-2,1))*dt;
Pos(i,2)=Pos(i-1,2)+1.5*(Pos(i-1,1)*(r-Pos(i-1,3))-Pos(i-1,2))*dt-0.5*(Pos(i-2,1)*(r-Pos(i-2,3))-Pos(i-2,2))*dt;
Pos(i,3)=Pos(i-1,3)+1.5*(Pos(i-1,1)*Pos(i-1,2)-b*Pos(i-1,3))*dt-0.5*(Pos(i-2,1)*Pos(i-2,2)-b*Pos(i-2,3))*dt;
end
head(Pos,10)