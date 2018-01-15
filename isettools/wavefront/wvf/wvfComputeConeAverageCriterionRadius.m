function [coneAvgCriterionRadius, wvfP] = wvfComputeConeAverageCriterionRadius(wvfP,defocus)
% Calculate the cone averaged criterion radius of the PSF
%
% Syntax:
%   [f, tmpWvfParams] = wvfComputeConeAverageCriterionRadius(wvfP,defocus)
%
% Description:
%    Calculate the cone averaged criterion radius of the PSF.
%    This calculation depends on a number of values specified in the
%    wavefront structure.
%       wvfGet(wvfP,'criterionFraction')
%       wvfGet(wvfP,'coneWeights');
%       wvfGet(wvfP,'- theThat is,
%    what is the radius that contains the criterion fraction of the PSF
%    mass, where that f
%
% Inputs:
%    wvfP            - Wavefront struct
%    defocus         - Defocus zernike coefficient, in diopters
%
% Outputs:
%    f            - Focus?
%    tmpWvfParams - wavefront struct with defocus stuck i..
%
% Optional key/value pairs
%     None.
%
% See also:
%

% History:
%    01/15/18  dhb  Started to bring this into the modern era


% Convert defocus from diopters to microns and set

wvfP = wvfSet(wvfP,'zcoeffs',defocus,'defocus');

% Compute the PSF with defocus set
wvfP = wvfComputeConePSF(wvfP);

% Get the criterion radius for each cone sensitivity and spectrum
% weighted PSF, and average these up, weighting the cones as
% specified.

nCones = size(wvfP.T_cones, 1);
criterionRadius = 0;
for j = 1:nCones
    criterionRadius(j) = psfFindCriterionRadius(...
        wvfP.conepsf(:, :, j), ...
        wvfP.criterionFraction);
    
    coneAvgCriterionRadius = coneAvgCriterionRadius + wvfP.coneWeights(j) * criterionRadius(j);
end

end
