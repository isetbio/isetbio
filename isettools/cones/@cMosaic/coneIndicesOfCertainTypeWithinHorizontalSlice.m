function [theConeIndices, theROI] = coneIndicesOfCertainTypeWithinHorizontalSlice(...
    obj, slicePositionDegs, sliceThicknessDegs, coneType)
% Method to return the indices of cones within a horizontal slice through the mosaic
%
% Syntax:
%   [theConeIndices, theROI] = coneIndicesOfCertainTypeWithinHorizontalSlice(...
%         obj, slicePositionDegs, sliceThicknessDegs, coneType)
%
% Description:
%    Method to return the indices of cones of a single specified cone type that 
%    lie within a horizontal slice through the mosaic
%
% Inputs:
%    obj                 - A @cMosaic object
%    slicePositionDegs   - Y-coord of slice
%    sliceThicknessDegs  - Y-thickness of slice 
%    coneType            - Type of cone. Choose from [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID]
%
% Outputs:                 
%    theConeIndices      - Indices of cones of specified type that lie 
%                          within the specified horizontal slice, sorted 
%                          in increasing x-coord eccentricity
%
% Optional key/value pairs:
%                        None

    % Define an ROI
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [obj.eccentricityDegs(1) slicePositionDegs], ...
            'width', obj.sizeDegs(1), ...
            'height', sliceThicknessDegs, ...
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