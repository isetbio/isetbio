function [visualizedConeIndices, mRGCIndices, theROI] = ...
    coneAndMRGCindicesAlongDesiredYposition(theMRGCMosaic, targetYdegs)

    % Define an ROI
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [theMRGCMosaic.eccentricityDegs(1) targetYdegs], ...
            'width', theMRGCMosaic.sizeDegs(1), ...
            'height', 0.1, ...
            'rotation', 0.0...
        ));

    visualizedConeIndices = theROI.indicesOfPointsInside(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs);
    mRGCIndices = theROI.indicesOfPointsInside(theMRGCMosaic.rgcRFpositionsDegs);

    theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
    theVisualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(mRGCIndices,1));

    [~,idx] = sort(theVisualizedConeXcoords, 'ascend');
    visualizedConeIndices = visualizedConeIndices(idx);

    % Exclude S-cones
    idx = find(theMRGCMosaic.inputConeMosaic.coneTypes(visualizedConeIndices) == cMosaic.SCONE_ID);
    [~, idx] = setdiff(visualizedConeIndices, visualizedConeIndices(idx));
    visualizedConeIndices = visualizedConeIndices(idx);
    
    [~,idx] = sort(theVisualizedMRGCXcoords, 'ascend');
    mRGCIndices = mRGCIndices(idx);
end