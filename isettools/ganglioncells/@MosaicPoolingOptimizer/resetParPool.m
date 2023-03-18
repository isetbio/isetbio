function  shutdownParPoolOnceCompleted = resetParPool(parPoolSize)
    if (~isempty(parPoolSize))
        poolobj = gcp('nocreate'); 
        if (~isempty(poolobj)) && ( poolobj.NumWorkers ~= parPoolSize)
           fprintf('Deleting previous parpool with %d workers to start a new one with %d workers instead\n', ...
               poolobj.NumWorkers, parPoolSize);
           delete(poolobj);
           shutdownParPoolOnceCompleted = true;
           parpool('local',parPoolSize);
        end
    end
end
