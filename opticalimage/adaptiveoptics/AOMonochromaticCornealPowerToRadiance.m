function radiance = AOMonochromaticCornealPowerToRadiance(wls,componentWavelengths,componentCornealPowerUW,pupilDiamMm,stimulusAreaDegs2)
% radiance = AOMonochromaticCornealPowerToRadiance(wls,componentWavelengths,componentCornealPowerUW,pupilDiamMm,stimulusAreaDegs2)
%
% Convert stimulus description in terms of monochromatic component
% wavelengths and their corneal power to an equivalent spectral radiance on
% wavelength sampling vector wls.
%
% Depends both on the area over which the power is spread in the pupil and
% the retinal area over which it is spread.
%
% Power is returned as Watts/[sr-m2-nm].

% Get pupil area
pupilAreaMm2 = pi*(pupilDiamMm/2)^2;
pupilAreaCm2 = pupilAreaMm2*10^-2;

% Background
nBgWavelengths = length(componentWavelengths);
radiance = zeros(length(wls),1);
for ww = 1:nBgWavelengths
    theWavelength = componentWavelengths(ww);
    theCornealIrradiance = componentCornealPowerUW(ww)/pupilAreaCm2;
    
    % UW is really UW/cm2 because the area of the detector is 1 cm2.  This
    % conversion gives us radiance in UW/[sr-cm2] for the narrowband laser
    % light.
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
