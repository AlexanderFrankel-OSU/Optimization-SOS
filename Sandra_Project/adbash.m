function [X] = adbash(X)
load("Sandra_Vars.mat","sig","b","r","dt")


    for i = 2
    X(i,1)=X(i-1,1)+sig*(X(i-1,2)-X(i-1,1))*dt;
    X(i,2)=X(i-1,2)+(X(i-1,1)*(r-X(i-1,3))-X(i-1,2))*dt;
    X(i,3)=X(i-1,3)+(X(i-1,1)*X(i-1,2)-b*X(i-1,3))*dt;
    end
    
    for i=3:length(X)
    X(i,1)=X(i-1,1)+1.5*(sig*(X(i-1,2)-X(i-1,1)))*dt-0.5*(sig*(X(i-2,2)-X(i-2,1)))*dt;
    X(i,2)=X(i-1,2)+1.5*(X(i-1,1)*(r-X(i-1,3))-X(i-1,2))*dt-0.5*(X(i-2,1)*(r-X(i-2,3))-X(i-2,2))*dt;
    X(i,3)=X(i-1,3)+1.5*(X(i-1,1)*X(i-1,2)-b*X(i-1,3))*dt-0.5*(X(i-2,1)*X(i-2,2)-b*X(i-2,3))*dt;
    end