function [coneAvgCriterionRadius, wvfP] = wvfComputeConeAverageCriterionRadius(wvfP,defocusDiopters)
% Calculate the cone averaged criterion radius of the PSF
%
% Syntax:
%   [f, tmpWvfParams] = wvfComputeConeAverageCriterionRadius(wvfP,defocusDiopters)
%
% Description:
%    Calculate the cone averaged criterion radius of the PSF.  This is the
%    radius that contains the criterion fraction of the PSF mass.
%
%    This calculation depends on a number of values specified in the
%    wavefront structure.
%       wvfGet(wvfP,'criterionFraction')
%       wvfGet(wvfP,'conePsfInfo');
%
% Inputs:
%    wvfP             - Input wavefront struct
%    defocusDiopters  - Defocus zernike coefficient, in diopters
%
% Outputs:
%    coneAvgCriterionRadius   - Cone averaged criterion radius
%    wvfPOut                  - Wavefront struct with passed defocus and
%                               psf computed.  Defocus is set via (not
%                               added to) wvfSet(wvfP,'calc observer
%                               accommodation'), not directly in the
%                               Zernike coefficients.
%
% Optional key/value pairs
%     None.
%
% See also:
%     wfvGet, conePsfInfoGet, wvfComputeOptimizedConePsf

% History:
%    01/15/18  dhb  Started to bring this into the modern era

% Examples:
%{

%}

% Set defocus via the observer accommodation field
wvfP = wvfSet(wvfP,'calc observer accommodation',defocusDiopters);

% Compute the PSF with defocus set
wvfP = wvfComputeConePSF(wvfP);

% Get the criterion radius for each cone sensitivity and spectrum
% weighted PSF, and average these up, weighting the cones as
% specified.
conePsfInfo = wvfGet(wvfP,'calc cone psf info');
T = conePsfInfoGet(conePsfInfo,'spectralSensitivities');
nCones = size(T, 1);
criterionRadius = 0;
for j = 1:nCones
    coneP
    criterionRadius(j) = psfFindCriterionRadius(...
        wvfP.conepsf(:, :, j), ...
        wvfP.criterionFraction);
    
    coneAvgCriterionRadius = coneAvgCriterionRadius + wvfP.coneWeights(j) * criterionRadius(j);
end

end
