% Method to resample the photocurrents using bicubic interpolation.
%
% NPC, ISETBIO Team 2016
% 

function resampledPhotocurrents = resample(photocurrents, originalTimeAxis, resampledTimeAxis)

    reshapeMatrix = false;
    
    if (ndims(photocurrents) == 3)
        % reshape for efficient computation
        [photocurrents , r, c] = RGB2XWFormat(photocurrents);
        reshapeMatrix = true;
    elseif (ndims(photocurrents) ~= 2)
        error('photocurrents must be either 2D or 3D');
    end
    
    % resample to the eye movement timebase, using bicubic intepolation
    resampledPhotocurrents = zeros(size(photocurrents, 1), numel(resampledTimeAxis));
    for coneIndex = 1:size(resampledPhotocurrents, 1)
        resampledPhotocurrents(coneIndex,:) = interp1(originalTimeAxis, squeeze(photocurrents(coneIndex,:)), resampledTimeAxis, 'pchip', 'extrap');
    end
    clear 'photocurrents'
    
    if (reshapeMatrix)
        resampledPhotocurrents = XW2RGBFormat(resampledPhotocurrents, r, c);
    end
end
