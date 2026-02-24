for i=3:500000
daPos(i,1)=daPos(i-1,1)+1.5*(sig*(daPos(i-1,2)-daPos(i-1,1))-ux*(daPos(i-1,1)-Pos(i-1,1)))*dt-0.5*(sig*(daPos(i-2,2)-daPos(i-2,1))-ux*(daPos(i-2,1)-Pos(i-2,1)))*dt;
daPos(i,2)=daPos(i-1,2)+1.5*((daPos(i-1,1)*(r-daPos(i-1,3))-daPos(i-1,2))-uy*(daPos(i-1,2)-Pos(i-1,2)))*dt-0.5*((daPos(i-2,1)*(r-daPos(i-2,3))-daPos(i-2,2))-uy*(daPos(i-2,2)-Pos(i-2,2)))*dt;
daPos(i,3)=daPos(i-1,3)+1.5*((daPos(i-1,1)*daPos(i-1,2)-b*daPos(i-1,3))-uz*(daPos(i-1,3)-Pos(i-1,3)))*dt-0.5*((daPos(i-2,1)*daPos(i-2,2)-b*daPos(i-2,3))-uz*(daPos(i-2,3)-Pos(i-2,3)))*dt;
end
head(daPos,10)