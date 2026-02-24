%This function takes in a given matrix and wipes everything but the initial
%values in the first row.
function [X] = initmatrix(X);
Width = size(X,2);
Height = size(X,1);

    for i = 2:Height
    X(i,:)=zeros(1,Width);
    end
end