function [theConeIndices, theROI] = ...
    coneIndicesOfCertainTypeWithinHorizontalSlice(obj, slicePositionDegs, sliceWidthDegs, coneType)

    % Define an ROI
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [obj.eccentricityDegs(1) slicePositionDegs], ...
            'width', obj.sizeDegs(1), ...
            'height', sliceWidthDegs, ...
            'rotation', 0.0...
        ));
    
    switch (coneType)
        case cMosaic.LCONE_ID
            idx = theROI.indicesOfPointsInside(obj.coneRFpositionsDegs(obj.lConeIndices,:));
            theConeIndices = obj.lConeIndices(idx);
        case cMosaic.MCONE_ID
            idx = theROI.indicesOfPointsInside(obj.coneRFpositionsDegs(obj.mConeIndices,:));
            theConeIndices = obj.mConeIndices(idx);
        case cMosaic.SCONE_ID
            idx = theROI.indicesOfPointsInside(obj.coneRFpositionsDegs(obj.sConeIndices,:));
            theConeIndices = obj.sConeIndices(idx);
    end

    if (isempty(theConeIndices))
        return;
    end

    theConeXcoords = squeeze(obj.coneRFpositionsDegs(theConeIndices,1));

    % Return cone indices sorted according to their x-coord
    [~,idx] = sort(theConeXcoords, 'ascend');
    theConeIndices = theConeIndices(idx);
end