function [M] = FrdEulCntrlTD(InitCons, Operator, Time, dt, Data, MU, MUprtr)
% Reminder: If the Operator is already just a matrix, you can turn this
% into a function handle and simplify everything greatly. That should cut
% down about half of this.
arguments
InitCons (:,:) double {mustBeNumeric}
Operator {mustbematrixorfunctionhandle}
Time double {mustBeNumeric}
dt double {mustBeNumeric}
Data (:,:) double = []
MU (:,:) double = []
MUprtr {mustbematrixorfunctionhandle} = []
end

N = ceil(Time/dt);
if size(InitCons,1)*N*8 > 1.41e+10
    M = NaN;
    return;
end

if isequal(InitCons,zeros(size(InitCons,1),1));
    display('Enter initial state')
   return;
end

M = [InitCons(:,1),zeros(size(InitCons,1),N)];

if isequal(Data,[]) % If no data assimilation, just run normally.
    switch num2str(isa(Operator,'double'))
        case '1'
            for col = 2:N+1
                    M(:,col) = M(:,col-1)+Operator*M(:,col-1)*dt;
            end
            return;
            
        case '0'
            if nargin(Operator)==1
                for col = 2:N+1
                    M(:,col) = M(:,col-1)+Operator(M(:,col-1))*dt;
                end
                return;
            elseif nargin(Operator)==2
                for col = 2:N+1
                    M(:,col) = M(:,col-1)+Operator(M(:,col-1),[])*dt;
                end
                return;
            elseif nargin(Operator)==3
                for col = 2:N+1
                    M(:,col) = M(:,col-1)+Operator(M(:,col-1),[],(col-1)*dt)*dt;
                end
                return;
            else
            display(['Operator should have 1,2, or 3 arguments if a function handle. Order is:' ...
                'Model, Data, Time.'])
            return;
            end
        otherwise
    end
end

signature = num2str([isa(Operator,'double'),isa(MUprtr,'double')]);

switch signature % Checks if the operator or mu operator are matrices or function handles
    case '1  1' % If both Operator and MUprtr are matrices

        for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator*M(:,col-1)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
        end
        return;
    
    case '1  0' % If the Operator is a matrix and the MUprtr is a function handle

        if ~(nargin(MUprtr)==3)
            display('Number of arguments for MUprtr should be 3, even if unused.')
            display('Order is: Model, Data, Time.')
            return;
        end

        for col = 2:N+1 
            M(:,col) = M(:,col-1)+(Operator*M(:,col-1)-MUprtr(M(:,col-1),Data(:,col-1),(col-1)*dt))*dt;
        end
        return;

    case '0  1' % If the Operator is a function handle and the MUprtr is a matrix
        if nargin(Operator)==1
            for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1))-MU*(M(:,col-1)-Data(:,col-1)))*dt;
            end
            return;
        elseif nargin(Operator)==2
            for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1),(col-1)*dt)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
            end
            return;
        elseif nargin(Operator)==3
            for col = 2:N+1 
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1),Data(:,col-1),(col-1)*dt)-MU*(M(:,col-1)-Data(:,col-1)))*dt;
            end
            return;
        end

    case '0  0' % If the Operator and MUprtr are both function handles
        % if or(~(nargin(Operator)==3),~(nargin(MUprtr)==3))
        %     display(['Problem: Please make sure that both Operator and MUprtr have 3 arguments if you want both ' ...
        %         'to be function handles, even if some are empty. You can add the ~ for a blank. The order ' ...
        %         'is: Model, Data, Time.'])
        %     return;
        % end
        if isequal(Data,[])
            display('To run data assimilation, input data.')
            return;
        else
            for col = 2:N+1
                M(:,col) = M(:,col-1)+(Operator(M(:,col-1),Data(:,col-1),(col-1)*dt)-MUprtr(M(:,col-1),Data(:,col-1),(col-1)*dt))*dt;
            end
        end

    otherwise
        display('Problem: None of the cases were triggered in Forward Euler.')
        return;
end

    
end


function mustbematrixorfunctionhandle(MATLABwhy)
        assert(isa(MATLABwhy,'double') || isa(MATLABwhy,'function_handle'),...
            'Must either be a matrix or function handle.')
end