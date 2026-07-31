function [TTC] = TTCfinder(Diff,dt)
    arguments
    Diff (:,:) double
    dt (1,1) double = 1 % If dt is not included, assume you want index.
    end

    ROWS = size(Diff,1);
    COLS = size(Diff,2);
    TTC = zeros(ROWS,1);

    try
        if anynan(Diff(:,end))
             notNAN = zeros(ROWS,1);
             startNAN = COLS; % If you find a NaN above, you can start looking much earlier!
             for row = 1:ROWS
                notNAN(row,1) = find(~isnan(Diff(row,1:startNAN)),1,"last");
                startNAN = min([notNAN(row,1),startNAN]);
             end
             notNANmin = min(notNAN);
             for row = 1:ROWS
                TTC(row,1) = (find(Diff(row,1:notNANmin),1,"last")+1)*dt;
             end

             if ~(max(TTC)<notNANmin*dt) % If the convergence doesn't happen until the NaN starts.
                 TTC = NaN;
                 return;
             end
             
        else
            for row = 1:ROWS
                TTC(row,1) = (find(Diff(row,:),1,"last")+1)*dt;
            end
        end

        if ~(max(TTC) < COLS*dt) % If the convergence doesn't happen in time.
            TTC = NaN;
        end
    catch
    display('Check TTCfinder definition.')
    end

    


end