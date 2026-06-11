function resampledPhotocurrents = resample(photocurrents, ...
    originalTimeAxis, resampledTimeAxis)
% Method to resample the photocurrents using bicubic interpolation.
%
% Syntax:
%   rexampledPhotocurrents = resample(photocurrents, originalTimeAxis, ...
%       resampledTimeAxis)
%
% Description:
%    Method to resample the photocurrents using bicubic interpolation.
%
% Inputs:
%    photoCurrents          - Original photocurrents
%    originalTimeAxis       - Original time axis
%    resampledTimeAxis      - The resampled time axis from which to pull
%                             the photocurrents.
%
% Outputs:
%    resampledPhotocurrents - The calculated resampled photocurrents on the
%                             resampled time axis.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  NPC  ISETBIO Team 2016
%    02/12/18  jnm  Formatting

    reshapeMatrix = false;
    
    if (ndims(photocurrents) == 3)
        % reshape for efficient computation
        [photocurrents , r, c] = RGB2XWFormat(photocurrents);
        reshapeMatrix = true;
    elseif (ndims(photocurrents) ~= 2)
        error('photocurrents must be either 2D or 3D');
    end
    
    % resample to the eye movement timebase, using bicubic intepolation
    resampledPhotocurrents = zeros(size(photocurrents, 1), ...
        numel(resampledTimeAxis));
    for coneIndex = 1:size(resampledPhotocurrents, 1)
        resampledPhotocurrents(coneIndex, :) = interp1(...
            originalTimeAxis, squeeze(photocurrents(coneIndex, :)), ...
            resampledTimeAxis, 'pchip', 'extrap');
    end
    clear 'photocurrents'
    
    if (reshapeMatrix)
        resampledPhotocurrents = XW2RGBFormat(...
            resampledPhotocurrents, r, c);
    end
end
