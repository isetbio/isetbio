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

    % Reset parpool
    shutdownParPoolOnceCompleted = MosaicPoolingOptimizer.resetParPool(parPoolSize);

    % Generate patterns
    H = zeros(nStim, size(X,1), size(X,2), 'single');

    if (parPoolSize ~= 0)
        parfor sIndex = 1:nStim
            fprintf('Hartley pattern %d of %d\n', sIndex, nStim);
            fx = lIndices(sIndex);
            fy = mIndices(sIndex);
            a = 2*pi*(fx*X + fy*Y);
            H(sIndex,:,:) = single(sin(a)+cos(a));
        end
    
        % Shutdown parpool
        if (shutdownParPoolOnceCompleted)
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
    else
        for sIndex = 1:nStim
            fprintf('Hartley pattern %d of %d\n', sIndex, nStim);
            fx = lIndices(sIndex);
            fy = mIndices(sIndex);
            a = 2*pi*(fx*X + fy*Y);
            H(sIndex,:,:) = single(sin(a)+cos(a));
        end
    end


    H = H / max(H(:));
    H = single(H);
end