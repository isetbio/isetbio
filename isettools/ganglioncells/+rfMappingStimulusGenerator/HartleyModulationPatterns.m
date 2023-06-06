function [H, lIndices, mIndices] = HartleyModulationPatterns(omega, stimSizeDegs, pixelSizeDegs, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizePatterns', false, @islogical);
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizePatterns = p.Results.visualizePatterns;

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
    [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool(parPoolSize);

    % Generate patterns
    H = zeros(nStim, size(X,1), size(X,2), 'single');

    if ((~isempty(parPoolSize)) && (parPoolSize>1)) || (isempty(parPoolSize))
        parfor sIndex = 1:nStim
            fprintf('Hartley pattern %d of %d (%d parallel workers)\n', sIndex, nStim, numWorkers);
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
            fprintf('Hartley pattern %d of %d (serially)\n', sIndex, nStim);
            fx = lIndices(sIndex);
            fy = mIndices(sIndex);
            a = 2*pi*(fx*X + fy*Y);
            H(sIndex,:,:) = single(sin(a)+cos(a));
            if (visualizePatterns)
                if (sIndex == 1)
                    hFig = figure(1);
                    set(hFig, 'Position', [10 10 2000 1000]);
                    ax1 = subplot(1,2,1);
                    ax2 = subplot(1,2,2);
                    colormap(gray(1024));
                    sfMax = 1/(2*pixelSizeDegs);
                    sfSupport = linspace(-sfMax, sfMax, pixelsNum);
                end
                imagesc(ax1,squeeze(H(sIndex,:,:)));
                axis(ax1, 'image');
                ftMag = fftshift(abs(fft2(squeeze(H(sIndex,:,:)))));
                imagesc(ax2,sfSupport, sfSupport, ftMag);
                axis(ax2, 'image');
                set(ax2, 'XLim', [-60 60], 'YLim', [-60 60], ...
                         'XTick', -60:10:60, 'YTick', -60:10:60);
                drawnow;
            end
        end
    end


    H = H / max(H(:));
    H = single(H);
end