function rfIndices = indicesOfRFsWithinROI(obj, geometryStruct)
% Method to return the indices of mRGC RFs within a geometryStruct
%
% Syntax:
%   rfIndices = obj.indicesOfRFsWithinROI(geometryStruct)
%
% Description:
%    return the indices of mRGC RFs within a geometryStruct appropriate for @regionOfInterest
%
% Inputs:
%    obj                 - An @mRGCmosaic object
%    geometryStruct      - A geometry struct appropriate for @regionOfInterest
%
% Outputs:                 Indices of cones within the region of interest

    theROI = regionOfInterest('geometryStruct', geometryStruct);

    switch (theROI.units)
        case 'microns'
            rfIndices = theROI.indicesOfPointsInside(obj.rgcRFpositionsMicrons);
        otherwise
            rfIndices = theROI.indicesOfPointsInside(obj.rgcRFpositionsDegs);
    end
end