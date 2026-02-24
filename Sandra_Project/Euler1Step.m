for i=2
Pos(i,1)=Pos(i-1,1)+sig*(Pos(i-1,2)-Pos(i-1,1))*dt;
Pos(i,2)=Pos(i-1,2)+(Pos(i-1,1)*(r-Pos(i-1,3))-Pos(i-1,2))*dt;
Pos(i,3)=Pos(i-1,3)+(Pos(i-1,1)*Pos(i-1,2)-b*Pos(i-1,3))*dt;
end
head(Pos,10)