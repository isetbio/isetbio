% Method to temporally resample an absorptions count sequence ensuring that the total
% number of photons are unchanged in the original and the resampled versions
function resampledAbsorptionsSequence = tResample(absorptionsSequence, pattern, originalTimeAxis, resampledTimeAxis)
 
    reshapeMatrix = false;
    
    if (ndims(absorptionsSequence) == 3)
        % reshape for efficient computation
        [absorptionsSequence , r, c] = RGB2XWFormat(absorptionsSequence);
        reshapeMatrix = true; 
    elseif (ndims(absorptionsSequence) ~= 2) && (ndims(absorptionsSequence) ~= 1)
        error('absorptionsCountSequence must be either 1D 2D or 3D');
    end
    
    if ((numel(originalTimeAxis) == 1) && (numel(resampledTimeAxis) == 1))
        resampledAbsorptionsSequence = absorptionsSequence;
        return;
    end
    
    % resample count sequence in the photocurrent timebase, i.e. upsample.
    resampledAbsorptionsSequence = zeros(size(absorptionsSequence,1), numel(resampledTimeAxis));
    for coneType = 2:4
        indices = find(pattern == coneType);
        for k = 1:numel(indices)
            coneIndex = indices(k);
            countsBefore = squeeze(absorptionsSequence(coneIndex,:));
            countsAfter  = interp1(originalTimeAxis, countsBefore, resampledTimeAxis, 'previous');
            resampledAbsorptionsSequence(coneIndex,:) = countsAfter / sum(countsAfter) * sum(countsBefore);
        end
    end
    
    if (reshapeMatrix)
        resampledAbsorptionsSequence = XW2RGBFormat(resampledAbsorptionsSequence, r, c);
    end
end