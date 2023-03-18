function [H, lIndices, mIndices] = HartleyModulationPatterns(omega, stimSizeDegs, pixelSizeDegs, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;

    pixelsNum  = round(stimSizeDegs / pixelSizeDegs);
    pixelIndex = 0:(pixelsNum-1);
    [X,Y] = meshgrid(pixelIndex/pixelsNum, pixelIndex/pixelsNum);

    nStim = (2*omega+1)^2;
   
    lIndices = zeros(nStim,1);
    mIndices = lIndices;
    sIndex = 0;
    for mIndex = 0:(2*omega)
        for lIndex = 0:(2*omega)
            sIndex = sIndex + 1;
            lIndices(sIndex) = lIndex-omega;
            mIndices(sIndex) = mIndex-omega;
        end
    end

    shutdownParPoolOnceCompleted = false;
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

    % Generate patterns
    H = zeros(nStim, size(X,1), size(X,2), 'single');
    parfor sIndex = 1:nStim
        fprintf('Hartley pattern %d of %d\n', sIndex, nStim);
        fx = lIndices(sIndex);
        fy = mIndices(sIndex);
        a = 2*pi*(fx*X + fy*Y);
        f = sin(a)+cos(a);
        H(sIndex,:,:) = single(f);
    end

    % Shutdown parpool
    if (shutdownParPoolOnceCompleted)
       delete(poolobj);
    end

    H = H / max(H(:));
    H = single(H);
end