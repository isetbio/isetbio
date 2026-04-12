%
% RGCMosaicConstructor.helper.simulateExperiment.updateBackgroundToAchieveEqualLandMconeActivation
%
function [backgroundChromaticity, backgroundLuminanceCdM2, LMSbefore, LMSafter] = ...
        updateBackgroundToAchieveEqualLandMconeActivation(...
            thePresentationDisplay, ...
            backgroundChromaticity, ...
            backgroundLuminanceCdM2 ...
        )

    fprintf(2,'Adjust backround chromaticity to achieve equal L and M cone excitations (2 deg Stockman only)\n')
    displayLinearRGBToXYZ = displayGet(thePresentationDisplay, 'rgb2xyz');
    displayLinearRGBToLMS = displayGet(thePresentationDisplay, 'rgb2lms');
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

