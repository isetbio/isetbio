function resampledAbsSequence = tResample(absSequence, pattern, originalTimes, resampledTimes)
%TRESAMPLE  Resample an absorption time sequence
%   resampledAbsorptions = tResample(absorptionsSequence, pattern, originalTimes, resampledTimes)
%
%   Resample an absorption sequence (x,y,t) to a new set of
%   sample times, (x,y,T). 
%
%   We ensure that the total number of photons are unchanged in the original
%   and the resampled versions
%
%   This is used to particularly with the osBioPhys conversion from absorptions to photocurrent.
%
%   Inputs:
%   [DHB COMMENT: PLEASE SAY WHAT THE INPUTS ARE.]
%
%   Outputs:
%   [DHB COMMENT: PLEASE SAY WHAT THE OUTPUTS ARE.]
%
%   Optional key/value pairs:
%   None.

% NC, ISETBIO Team, 2016

%% Check that input is actually a time sequence.
if (ndims(absSequence) == 3)
    % reshape for efficient computation
    [absSequence , r, c] = RGB2XWFormat(absSequence);
else
    error('absorptions must vary over time (3D)');
end

% Should never be in this case because we are checking above that the absorptions vary over time.
% But what the heck, handle it anyway.
if ((numel(originalTimes) == 1) && (numel(resampledTimes) == 1))
    resampledAbsSequence = absSequence;
    return;
end

%% Resample count sequence in the photocurrent timebase, i.e. upsample.
resampledAbsSequence = zeros(size(absSequence,1), numel(resampledTimes));
for coneType = 2:4
    indices = find(pattern == coneType);
    for k = 1:numel(indices)
        coneIndex = indices(k);
        countsBefore = squeeze(absSequence(coneIndex,:));
        countsAfter  = interp1(originalTimes, countsBefore, resampledTimes, 'previous');
        resampledAbsSequence(coneIndex,:) = countsAfter / sum(countsAfter) * sum(countsBefore);
    end
end
resampledAbsSequence = XW2RGBFormat(resampledAbsSequence, r, c);

end