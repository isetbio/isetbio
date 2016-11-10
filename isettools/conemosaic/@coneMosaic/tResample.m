function resampledAbsorptionsSequence = tResample(absorptionsSequence, pattern, originalTimeAxis, resampledTimeAxis)
% Resample the absorption time series
%
% Method to temporally resample an absorptions count sequence ensuring that
% the total number of photons are unchanged in the original and the
% resampled versions
%
% NC, ISETBIO Team, 2016


if (ndims(absorptionsSequence) == 3)
    % reshape for efficient computation
    [absorptionsSequence , r, c] = RGB2XWFormat(absorptionsSequence);
else
    error('absorptions must vary over time (3D)');
end

% Should never be in this case because we are checking that the absorptions
% vary over time
if ((numel(originalTimeAxis) == 1) && (numel(resampledTimeAxis) == 1))
    resampledAbsorptionsSequence = absorptionsSequence;
    return;
end

% Resample count sequence in the photocurrent timebase, i.e. upsample.
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

resampledAbsorptionsSequence = XW2RGBFormat(resampledAbsorptionsSequence, r, c);

end