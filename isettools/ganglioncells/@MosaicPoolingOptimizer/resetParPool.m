function  [shutdownParPoolOnceCompleted, numWorkers] = resetParPool(parPoolSize)
    shutdownParPoolOnceCompleted = [];
    if ((~isempty(parPoolSize)) && (parPoolSize>1)) || (isempty(parPoolSize))
        poolobj = gcp('nocreate'); 
        if (~isempty(poolobj))
            numWorkers = poolobj.NumWorkers;
            if (numWorkers ~= parPoolSize)
               fprintf('Deleting previous parpool with %d workers to start a new one with %d workers instead\n', ...
                   poolobj.NumWorkers, parPoolSize);
               delete(poolobj);
               shutdownParPoolOnceCompleted = true;
               parpool('local',parPoolSize);
            end
        else
            parpool;
        end
    end

    poolobj = gcp('nocreate');
    if (isempty(poolobj))
        if ((~isempty(parPoolSize)) && (parPoolSize>1))
            parpool('local',parPoolSize);
        else
            parpool;
        end
    end
    poolobj = gcp('nocreate');
    numWorkers = poolobj.NumWorkers;
end
