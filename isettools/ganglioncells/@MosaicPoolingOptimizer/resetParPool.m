function  [shutdownParPoolOnceCompleted, numWorkers] = resetParPool(parPoolSize)
    shutdownParPoolOnceCompleted = false;

    if ((~isempty(parPoolSize)) && (parPoolSize>1))
        poolobj = gcp('nocreate'); 
        if (~isempty(poolobj))
            numWorkers = poolobj.NumWorkers;
            if (numWorkers ~= parPoolSize)
               delete(poolobj);
               fprintf('Deleting previous parpool with %d workers to start a new one with %d workers instead\n', ...
                   poolobj.NumWorkers, parPoolSize);
               pause(1.0);
               parpool('local',parPoolSize);
            end
        else
            parpool;
        end
    end

    poolobj = gcp('nocreate');
    numWorkers = poolobj.NumWorkers;
    fprintf('Number of parallel workers: %d\n', numWorkers);
end
