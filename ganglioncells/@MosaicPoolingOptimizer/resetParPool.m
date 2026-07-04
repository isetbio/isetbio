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
    if isempty(poolobj)
        numWorkers = 0;
        fprintf('No parallel pool is running.\n');
        return;
    end

    if ~isprop(poolobj, 'NumWorkers')
        numWorkers = 0;
        fprintf('The active parallel pool does not report a worker count.\n');
        return;
    end

    numWorkers = poolobj.NumWorkers;
    fprintf('Number of parallel workers: %d\n', numWorkers);
end
