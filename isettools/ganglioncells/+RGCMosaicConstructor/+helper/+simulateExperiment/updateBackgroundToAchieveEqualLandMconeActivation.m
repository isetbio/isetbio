%
% RGCMosaicConstructor.helper.simulateExperiment.updateBackgroundToAchieveEqualLandMconeActivation
%
function [backgroundChromaticity, backgroundLuminanceCdM2, LMSbefore, LMSafter] = ...
        updateBackgroundToAchieveEqualLandMconeActivation(...
            thePresentationDisplay, ...
            backgroundChromaticity, ...
            backgroundLuminanceCdM2, ...
            varargin)

    p = inputParser;
    p.addParameter('coneFundamentalsToEmploy', []);
    p.parse(varargin{:});

    coneFundamentalsToEmploy = p.Results.coneFundamentalsToEmploy;

    if (isempty(coneFundamentalsToEmploy))
        fprintf(2,'Adjust backround chromaticity to achieve equal L and M cone excitations (2 deg SS cone fundamentals)\n');
        displayLinearRGBToLMS = displayGet(thePresentationDisplay, 'rgb2lms');
    else
        fprintf(2,'Adjust backround chromaticity to achieve equal L and M cone excitations (custom cone fundamentals)\n');

        customConeFundamentals = coneFundamentalsToEmploy;
        assert(isfield(customConeFundamentals, 'wavelengthSupport'), ...
                'customConeFundamentals does not contain wavelength support info');
        assert(isfield(customConeFundamentals, 'quantalExcitationSpectra'), ...
                'customConeFundamentals does not contain quantalExcitationSpectra info');
        assert(size(customConeFundamentals.quantalExcitationSpectra,2) == 3, ...
                'customConeFundamentals.spd is not an Nx3 matrix');
        assert(size(customConeFundamentals.quantalExcitationSpectra,1) == numel(customConeFundamentals.wavelengthSupport), ...
                'customConeFundamentals.spf does not have the same dimensionality as customConeFundamentals.wavelengthSupport');
    
        % Match customConeFundamentals
        displayWavelengths = displayGet(thePresentationDisplay, 'wave');
        if (~isequal(displayWavelengths, customConeFundamentals.wavelengthSupport)) || ...
           ((isequal(displayWavelengths, customConeFundamentals.wavelengthSupport))&&(~all(displayWavelengths==customConeFundamentals.wavelengthSupport)))
            % Resample customConeFundamentals.spd to wavelength support matching the display
            resampledCustomConeFundamantals = displayWavelengths*0;
            for iChannel = 1:size(customConeFundamentals.quantalExcitationSpectra,2)
                resampledCustomConeFundamantals(:,iChannel) = interp1(...
                    customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,iChannel), ...
                    displayWavelengths, 'linear','extrap');
            end
            customConeFundamentals.quantalExcitationSpectra = resampledCustomConeFundamantals;
            customConeFundamentals.wavelengthSupport = displayWavelengths;
        end

        coneFundamentals = customConeFundamentals.quantalExcitationSpectra/max(customConeFundamentals.quantalExcitationSpectra(:));

        % Compute the displayRGCtoLMS matrix
        displayLinearRGBToLMS = (coneFundamentals' * displayGet(thePresentationDisplay, 'spd', displayWavelengths))';
    end


    displayLinearRGBToXYZ = displayGet(thePresentationDisplay, 'rgb2xyz');
    
    displayLMSToLinearRGB = inv(displayLinearRGBToLMS);
    displayXYZtoLinearRGB = inv(displayLinearRGBToXYZ);
    
    XYZ = xyYToXYZ([backgroundChromaticity(1) backgroundChromaticity(2)  backgroundLuminanceCdM2]');
    RGB = imageLinearTransform(XYZ', displayXYZtoLinearRGB);
    LMSbefore = imageLinearTransform(RGB, displayLinearRGBToLMS);
    
    backgroundLMSconeExcitations = LMSbefore;
    backgroundLMSconeExcitations(1:2) = 0.5*sum(LMSbefore(1:2));
    
    backgroundRGB = imageLinearTransform(backgroundLMSconeExcitations, displayLMSToLinearRGB);
    backgroundXYZ = imageLinearTransform(backgroundRGB, displayLinearRGBToXYZ);
    
    xyY = XYZToxyY(backgroundXYZ(:));
    LMSafter = imageLinearTransform(backgroundRGB, displayLinearRGBToLMS);
    
    backgroundChromaticity(1) = xyY(1);
    backgroundChromaticity(2) = xyY(2);
    backgroundLuminanceCdM2 = xyY(3);
end

