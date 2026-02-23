
%Argument 1 is the 3-column matrix to perform Adams-Bashform on. Argument 2
%is the data to use for data assimilation. They should have the same size.

function [X] = loradbashDA(X,Y)
load("Sandra_Vars.mat","sig","b","r","ux","uy","uz", "dt", "Pos")
MatWidth = size(X,2);

if (X(1,:) == zeros(1,MatWidth))
    print("Please enter initial conditions.")

elseif size(X)~=size(Y)
    print("Size of given matrix and data are different.")

else
    for i = 2
    X(i,1)=X(i-1,1)+(sig*(X(i-1,2)-X(i-1,1))-ux*(X(i-1,1)-Y(i-1,1)))*dt;
    X(i,2)=X(i-1,2)+((X(i-1,1)*(r-X(i-1,3))-X(i-1,2))-uy*(X(i-1,2)-Y(i-1,2)))*dt;
    X(i,3)=X(i-1,3)+(X(i-1,1)*X(i-1,2)-b*X(i-1,3)-uz*(X(i-1,3)-Y(i-1,3)))*dt;
    end
    
    for i=3:length(X)
    X(i,1)=X(i-1,1)+1.5*(sig*(X(i-1,2)-X(i-1,1))-ux*(X(i-1,1)-Pos(i-1,1)))*dt-0.5*(sig*(X(i-2,2)-X(i-2,1))-ux*(X(i-2,1)-Y(i-2,1)))*dt;
    X(i,2)=X(i-1,2)+1.5*((X(i-1,1)*(r-X(i-1,3))-X(i-1,2))-uy*(X(i-1,2)-Pos(i-1,2)))*dt-0.5*((X(i-2,1)*(r-X(i-2,3))-X(i-2,2))-uy*(X(i-2,2)-Y(i-2,2)))*dt;
    X(i,3)=X(i-1,3)+1.5*((X(i-1,1)*X(i-1,2)-b*X(i-1,3))-uz*(X(i-1,3)-Pos(i-1,3)))*dt-0.5*((X(i-2,1)*X(i-2,2)-b*X(i-2,3))-uz*(X(i-2,3)-Y(i-2,3)))*dt;
    end
    
end



