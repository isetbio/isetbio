function [coneSpacingDegs, coneSpacingMMs, coneDensityDegs2, coneDensityMMs2] = ...
    coneSpacingAtRetinalPositions(whichEye, conePosMicrons)
% Compute cone spacing and densities for arbitrary retinal positions in either eye.

    % Check size of input
    assert(size(conePosMicrons,2) == 2, sprintf('conePosMicrons must be an N x 2 matrix of (x,y) positions'));
    
    % Compute ecc radii and angles
    eccRadiiMicrons  = sqrt(sum(conePosMicrons.^2,2));
    eccRetinalAngles = atan2d(conePosMicrons(:,2),conePosMicrons(:,1));
    
    % Reshape to [1 X N]
    eccRadiiMicrons = reshape(eccRadiiMicrons, [1 numel(eccRadiiMicrons)]);
    eccRetinalAngles = reshape(eccRetinalAngles, [1 numel(eccRetinalAngles)]);
    
    % Convert retinal angles at whichEye to their equivalent right eye visual field angles
    rightEyeEccVisualAngles = RGCmodels.Watson.convert.retinalAnglesAtSpecificEyeToAnglesInRightVisualField(eccRetinalAngles, whichEye);
    
    % Convert ecc radii in microns to degs
    eccRadiiDegs = RGCmodels.Watson.convert.rhoMMsToDegs(eccRadiiMicrons * 1e-3);

    % Compute cone densities along right eye visual meridians for all eccRadii
    indexedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    meridianConeDensitiesMMs2 = zeros(numel(indexedMeridians), numel(eccRadiiDegs));
    for meridianIndex = 1:numel(indexedMeridians)
        rightEyeVisualFieldMeridianName = indexedMeridians{meridianIndex};
        fprintf('meridian %d for %s: %s\n', meridianIndex, whichEye, rightEyeVisualFieldMeridianName);
        [~,~, ~,meridianConeDensitiesMMs2(meridianIndex,:)] = ...
            RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccRadiiDegs, rightEyeVisualFieldMeridianName);
    end
    
    % Interpolate densities from the meridian densities
    interpolationMethod = 'makima'; % 'linear'; % 'spline'; % 'linear'
    coneDensityMMs2 = RGCmodels.Watson.compute.radiallyInterpolated2DMapFromMeridianValues(...
        meridianConeDensitiesMMs2, rightEyeEccVisualAngles, interpolationMethod);
    
    % Convert density per mm^2 to density per deg^2
    coneDensityDegs2 = RGCmodels.Watson.convert.densityMMs2ToDegs2(coneDensityMMs2, eccRadiiDegs);
    
    % Convert density to spacing in degs
    coneSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(coneDensityDegs2); 
    
    % Convert density to spacing in MMs
    coneSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(coneDensityMMs2); 
end

