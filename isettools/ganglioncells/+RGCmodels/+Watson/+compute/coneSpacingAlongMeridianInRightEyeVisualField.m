function [coneSpacingDegs, coneSpacingMMs, coneDensityDegs2, coneDensityMMs2] = coneSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName)
% Compute cone spacing&density for desired cone positions (specified in visual degrees) along a meridian on the right eye visual field

    % Ensure all eccentricities are positive
    assert(all(eccDegs>=0), sprintf('Eccentrity vector must not include any negative values'));
    
    % Determine the isetbio angle (in the right eye) that corresponds to the rightEyeVisualFieldMeridianName
    isetbioAngle = RGCmodels.Watson.convert.rightEyeVisualFieldMeridianToISETBioAngleInRightEye(rightEyeVisualFieldMeridianName);
    
    % Convert eccDegs to eccMMs so that we don't go through the baked-in
    % ecc-independent 0.3 mm / deg conversion in coneDensityReadData()
    eccMM = RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    
    % Read raw Curcio dendity data
    [~, ~, coneDensityMMs2] = coneSizeReadData('eccentricity', eccMM, ...
                                        'angle', isetbioAngle*ones(1,numel(eccMM)), ...
                                        'eccentricityUnits', 'mm', ...
                                        'angleUnits','deg', ...
                                        'whichEye', 'right', ...
                                        'useParfor', true);
               

    % Apply connection for ISETBio foveal cone density being a bit higher (18,800 cones/deg2 vs 14,804 in Watson)                                
    coneDensityMMs2 = RGCmodels.Watson.compute.correctionForISETBioFovealConeDensity(coneDensityMMs2, eccDegs);

    % Cone density to spacing in MMs
    coneSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(coneDensityMMs2);
    
    % Convert density per mm^2 to density per deg^2
    coneDensityDegs2 = RGCmodels.Watson.convert.densityMMs2ToDegs2(coneDensityMMs2, eccDegs);
    
    % Convert density to spacing in degs
    coneSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(coneDensityDegs2); 
end

