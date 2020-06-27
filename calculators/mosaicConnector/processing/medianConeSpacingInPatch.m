function [coneSpacingMicrons, coneSpacingMaxMicrons] = medianConeSpacingInPatch(whichEye, patchEccMicrons, patchSizeMicrons)

    w = WatsonRGCModel('generateAllFigures', false);
    posUnits = 'mm'; densityUnits = 'mm^2';
    
    % Cone spacing at the center of the patch
    coneSpacingMicronsAtPatchCenter = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
        patchEccMicrons*1e-3, whichEye, posUnits, densityUnits, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
    % Generate grid of nodes on which to evaluare cone spacing
    xx = 0:coneSpacingMicronsAtPatchCenter:(0.5*patchSizeMicrons(1));
    xx = [-fliplr(xx) xx(2:end)];
    x = patchEccMicrons(1) + xx;
    yy = 0:coneSpacingMicronsAtPatchCenter:(0.5*patchSizeMicrons(2));
    yy = [-fliplr(yy) yy(2:end)];
    y = patchEccMicrons(2) + yy;
    [X,Y] = meshgrid(x,y);
    X = X(:); Y = Y(:);

    % Determine cone spacings at all nodes of the grid
    coneSpacingMicronsList = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
            [X Y]*1e-3, whichEye, posUnits, densityUnits, ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
    % Determine the median cone spacing with the patch
    coneSpacingMicrons = median(coneSpacingMicronsList);
    coneSpacingMaxMicrons = max(coneSpacingMicronsList);
    fprintf('Cone spacing (patch center): %2.1f microns, median over patch: %2.1f\n', ...
        coneSpacingMicronsAtPatchCenter, coneSpacingMicrons);
end

