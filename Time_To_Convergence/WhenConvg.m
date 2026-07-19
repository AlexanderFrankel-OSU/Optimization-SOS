function [IndCnvg] = WhenConvg(A,B,CnvgPer)
% This script just finds where two matrices agree for a period of CnvgPer.
% Runs faster with larger CnvgPer.
% Checks and balances (:
arguments
    A (:,:) double
    B (:,:) double
    CnvgPer (1,1) double = 1
    
end

if ne(size(A),size(B))
    display('Error. Different-sized arrays.')
    return;
end

if CnvgPer < 1
    display('Positive natural number convergence period, please.')
    return;
end

if size(A,1) > size(A,2)
A = transpose(A); % Bad temporary solution.
B = transpose(B);
end

if ne(CnvgPer,1)
    CnvgPer = ceil(CnvgPer);
end

EvoLen = max(size(A));

if ne(CnvgPer,1)
    StepInts = ceil(EvoLen/abs(CnvgPer));
else
    StepInts = 1;
end

if CnvgPer == 1
    for ChckInt = 1:EvoLen;
        if isequal(A(:,ChckInt),B(:,ChckInt))
            IndCnvg = ChckInt;
            return;
        end
    end

else 
    % I.e., if we're not just running through every single index, 
    % then do it with ChckInt as an actual index matrix.
    ChckInt = 1:CnvgPer;
    while ChckInt(1) < (StepInts-1)*CnvgPer
        
        if isequal(A(:,ChckInt),B(:,ChckInt))
        IndCnvg = ChckInt(1);
        return;
        % At this point, if you want higher precision as to where the
        % convergence starts, we need to backtrack one interval and then 
        % check a bunch of other stuff.. It would be best to avoid this..

        end
        ChckInt = ChckInt+CnvgPer;
    end

    ChckInt = ChckInt(1):EvoLen;

    if isequal(A(:,ChckInt),B(:,ChckInt))
        IndCnvg = ChckInt(1);
        return;
    end

end

display('Does not converge according to the stepsize and tolerance.')
IndCnvg = NaN;
end
