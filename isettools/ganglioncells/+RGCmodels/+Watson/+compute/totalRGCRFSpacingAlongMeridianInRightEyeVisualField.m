function [rgcSpacingDegs, rgcSpacingMMs, rgcDensityDegs2, rgcDensityMMs2] = totalRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName)
% Compute total RGC receptive field spacing&density for desired positions (specified in visual degrees) along a meridian on the right eye visual field

    % Ensure all eccentricities are positive
    assert(all(eccDegs>=0), sprintf('Eccentrity vector must not include any negative values'));
    
    % Compute peak cone density in Degs2 and MMs2
    [~,~, coneDensityDegs2Peak, coneDensityMMs2Peak] = ...
        RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(...
        0.0,rightEyeVisualFieldMeridianName);

    % The foveal density of midget RGC RFs is twice the cone density
    % because in the fovea each cone connects to exactly 2 midget RGCs (one
    % ON and one OFF).
    midgetRGCRFDensityDegs2Peak = 2 * coneDensityDegs2Peak;
    
    % Retrieve percentage of total ganglion cells that are midget at 0 eccentricity
    percentMidgetsAtZeroEccentricity = RGCmodels.Watson.constants.f0;
    
    % Foveal density of total number of RGCs 
    totalRGCRFDensityDegs2Peak = 1.0/percentMidgetsAtZeroEccentricity * midgetRGCRFDensityDegs2Peak;
    
    % Retrieve params for the eccentricity variation along the rightEyeVisualFieldMeridianName
    mParams = RGCmodels.Watson.constants.meridianParamsTable(rightEyeVisualFieldMeridianName);
    
    % Compute normalized variation of RGC density with eccentricity for the rightEyeVisualFieldMeridianName
    eccVariation = RGCmodels.Watson.constants.totalRGCRFDensityEccVariation(...
        eccDegs, [mParams.a_k, mParams.r_2k, mParams.r_ek]);
    
    % Total RGC RF density per deg2, along the requested meridian
    rgcDensityDegs2 = totalRGCRFDensityDegs2Peak * eccVariation;
    
    % Convert density to spacing in degs
    rgcSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityDegs2); 
    
    % Repeat for density per MMs2
    midgetRGCRFDensityMMs2Peak = 2 * coneDensityMMs2Peak;
    totalRGCRFDensityMMs2Peak = 1.0/percentMidgetsAtZeroEccentricity * midgetRGCRFDensityMMs2Peak;
    
    % Total RGC RF density per MMs2, along the requested meridian
    rgcDensityMMs2 = totalRGCRFDensityMMs2Peak * eccVariation;
    
    % Convert density to spacing in MMs
    rgcSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityMMs2); 
end


