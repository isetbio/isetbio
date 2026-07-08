function [rgcSpacingDegs, rgcSpacingMMs, rgcDensityDegs2, rgcDensityMMs2] = midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName)
% Compute midget RGC receptive field spacing&density for desired positions (specified in visual degrees) along a meridian on the right eye visual field
% Note: This is for the combination of ON and OFF subtypes

    % Compute total RGC receptive field spacing&density
    [~, ~, totalrgcDensityDegs2, totalrgcDensityMMs2] = ...
        RGCmodels.Watson.compute.totalRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName);
    
    % Multiply with variation of midget RGC/totalRGC ratio with eccentricity
    rgcDensityDegs2 = RGCmodels.Watson.constants.midgetRGCFractionEccVariation(eccDegs) .* totalrgcDensityDegs2;
    rgcDensityMMs2  = RGCmodels.Watson.constants.midgetRGCFractionEccVariation(eccDegs) .* totalrgcDensityMMs2;
    
    % Convert density to spacing in degs
    rgcSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityDegs2); 
    rgcSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityMMs2); 
end


