function radiance = AOMonochromaticCornealPowerToRadiance(wls,componentWavelengths,componentCornealPowerUW,pupilDiamMm,stimulusAreaDegs2)
% Convert power at cornea to equivalent radiance
%
% Synopsis:
%   radiance = AOMonochromaticCornealPowerToRadiance(wls,componentWavelengths,componentCornealPowerUW,pupilDiamMm,stimulusAreaDegs2)
%
% Description:
%   Convert stimulus description in terms of narrowband stimulus
%   components, each centered on a specified wavelength, and their corneal
%   power to an equivalent spectral radiance on wavelength sampling vector
%   wls.
%
%   The meausurement of corneal power should be the total power that passes
%   through the pupil. This routine converts that to corneal irradiance by
%   dividing the passed number by the pupil area.
% 
%   The equivalent radiance is specific to the specified pupil and retinal
%   area, and is calculated to produce the same retinal illuminance as the
%   measured stimuli, given these two numbers.
% 
%   Power is passed as uW at each monochromatic wavelength specified, and
%   is returned as Watts/[sr-m2-nm].
%
% Inputs:
%   wls                           - A set of evenly space wavelengths
%   componentWavelengths          - The wavelengths where monochromatic power was measured. 
%                                   These must be a subset of the passed wls.
%   componentCornealPowerUW       - The corresponding powers in uW.
%   pupilDiamMm                   - Pupil diameter in mm.
%   stimulusAreaDegs2             - Stimulus area on retina in degrees squared.
%
% Outputs:
%   radiance                      - Radiance in Watts/[sr-m2-nm] on the passed wls support.
%
% Optional key/value pairs:
%   None.

% History:
%    12/13/20  dhb  Moved into ISETBio. 
%              dhb  ISETBio style comments.

% Get pupil area
pupilAreaMm2 = pi*(pupilDiamMm/2)^2;
pupilAreaCm2 = pupilAreaMm2*10^-2;

% Background
nBgWavelengths = length(componentWavelengths);
radiance = zeros(length(wls),1);
for ww = 1:nBgWavelengths
    theWavelength = componentWavelengths(ww);
    
    % Convert to corneal irradiance assuming that all the power passes
    % through the pupil.
    theCornealIrradiance = componentCornealPowerUW(ww)/pupilAreaCm2;
    
    % Convert to get radiance in UW/[sr-cm2] for the narrowband light at
    % this wavelength.
    radianceRaw = CornIrradianceAndDegrees2ToRadiance(theCornealIrradiance,stimulusAreaDegs2);
    
    % Convert to Watts/[sr-m2-nm] where we take the wavelength sampling
    % into account so that in the end the calculation of cone responses
    % will come out correctly.
    index = find(theWavelength == wls);
    if (length(index) ~= 1)
        error('Something funky about wls');
    end
    radiance(index) = (10^4)*(10^-6)*radianceRaw/(wls(2)-wls(1));
end
        
        
end
