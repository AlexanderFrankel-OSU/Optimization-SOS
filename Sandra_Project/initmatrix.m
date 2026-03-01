%This function takes in a given matrix and wipes everything but the initial
%values in the first row.
function [X] = initmatrix(X);
Width = size(X,2);
Height = size(X,1);

X=[X(1,:);zeros(Height-1,Width)];

end