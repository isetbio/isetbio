function coneIndices = indicesOfConesWithinROI(obj, geometryStruct)
% Method to return the indices of cones within a geometryStruct
%
% Syntax:
%   coneIndices = obj.indicesOfConesWithinROI(geometryStruct)
%
% Description:
%    return the indices of cones within a geometryStruct appropriate for @regionOfInterest
%
% Inputs:
%    obj                 - A @cMosaic object
%    geometryStruct      - A geometry struct appropriate for @regionOfInterest
%
% Outputs:                 Indices of cones within the region of interest

    theROI = regionOfInterest('geometryStruct', geometryStruct);

    switch (theROI.units)
        case 'microns'
            coneIndices = theROI.indicesOfPointsInside(obj.coneRFpositionsMicrons);
        otherwise
            coneIndices = theROI.indicesOfPointsInside(obj.coneRFpositionsDegs);
    end
end
