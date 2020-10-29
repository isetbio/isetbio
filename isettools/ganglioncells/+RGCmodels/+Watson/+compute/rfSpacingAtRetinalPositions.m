function [rfSpacingDegs, rfSpacingMMs, rfDensityDegs2, rfDensityMMs2] = ...
    rfSpacingAtRetinalPositions(whichEye, rfPosMicrons, neuronType)
% Compute rf spacing for specific neuron types and densities for arbitrary retinal positions in either eye.

    % Check size of input
    assert(size(rfPosMicrons,2) == 2, sprintf('conePosMicrons must be an N x 2 matrix of (x,y) positions'));
    
    % Compute ecc radii and angles
    eccRadiiMicrons  = sqrt(sum(rfPosMicrons.^2,2));
    eccRetinalAngles = atan2d(rfPosMicrons(:,2), rfPosMicrons(:,1));
    
    % Reshape to [1 X N]
    eccRadiiMicrons = reshape(eccRadiiMicrons, [1 numel(eccRadiiMicrons)]);
    eccRetinalAngles = reshape(eccRetinalAngles, [1 numel(eccRetinalAngles)]);
    
    % Convert retinal angles at whichEye to their equivalent right eye visual field angles
    rightEyeEccVisualAngles = RGCmodels.Watson.convert.retinalAnglesAtSpecificEyeToAnglesInRightVisualField(eccRetinalAngles, whichEye);
    
    % Convert ecc radii in microns to degs
    eccRadiiDegs = RGCmodels.Watson.convert.rhoMMsToDegs(eccRadiiMicrons * 1e-3);

    % Compute cone densities along right eye visual meridians for all eccRadii
    indexedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    meridianDensitiesMMs2 = zeros(numel(indexedMeridians), numel(eccRadiiDegs));
    
    for meridianIndex = 1:numel(indexedMeridians)
        rightEyeVisualFieldMeridianName = indexedMeridians{meridianIndex};
        
        switch (neuronType)
            case 'cones'
                [~,~, ~,meridianDensitiesMMs2(meridianIndex,:)] = ...
                    RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccRadiiDegs, rightEyeVisualFieldMeridianName);
            case 'all ganglion cells'
                [~,~, ~,meridianDensitiesMMs2(meridianIndex,:)] = ...
                    RGCmodels.Watson.compute.totalRGCRFSpacingAlongMeridianInRightEyeVisualField(eccRadiiDegs, rightEyeVisualFieldMeridianName);
            case 'midget ganglion cells'
                [~,~, ~,meridianDensitiesMMs2(meridianIndex,:)] = ...
                    RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(eccRadiiDegs, rightEyeVisualFieldMeridianName);
            otherwise
                error('Unknown neuronType: ''%s''.', neuronType);
        end
            
    end
    
    % Interpolate densities from the meridian densities
    interpolationMethod = 'makima'; % 'linear'; % 'spline'; % 'linear'
    rfDensityMMs2 = RGCmodels.Watson.compute.radiallyInterpolated2DMapFromMeridianValues(...
        meridianDensitiesMMs2, rightEyeEccVisualAngles, interpolationMethod);
    
    % Convert density per mm^2 to density per deg^2
    rfDensityDegs2 = RGCmodels.Watson.convert.densityMMs2ToDegs2(rfDensityMMs2, eccRadiiDegs);
    
    % Convert density to spacing in degs
    rfSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rfDensityDegs2); 
    
    % Convert density to spacing in MMs
    rfSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rfDensityMMs2); 
end

