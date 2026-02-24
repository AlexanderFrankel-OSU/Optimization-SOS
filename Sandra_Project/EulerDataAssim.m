for i=2
daPos(i,1)=daPos(i-1,1)+(sig*(daPos(i-1,2)-daPos(i-1,1))-ux*(daPos(i-1,1)-Pos(i-1,1)))*dt;
daPos(i,2)=daPos(i-1,2)+((daPos(i-1,1)*(r-daPos(i-1,3))-daPos(i-1,2))-uy*(daPos(i-1,2)-Pos(i-1,2)))*dt;
daPos(i,3)=daPos(i-1,3)+(daPos(i-1,1)*daPos(i-1,2)-b*daPos(i-1,3)-uz*(daPos(i-1,3)-Pos(i-1,3)))*dt;
end
head(daPos,10)