function coneIndices = indicesOfConesWithinROI(obj, roi)
% Method to return the indices of cones within a region of interest
%
% Syntax:
%   coneIndices = obj.indicesOfConesWithinROI(roi)
%
% Description:
%    return the indices of cones within a region of interest. The input
%    argument roi is a struct which specified either a rectangle or an
%    ellipse.
%
% Inputs:
%    obj                 - A @cMosaic object
%    roi                 - A struct specifiying either a rectangle or an
%                          ellipse
%
% Outputs:                 Indices of cones within the region of interest

    opticDiskROI = regionOfInterest('geometryStruct', roi);

    switch (opticDiskROI.units)
        case 'microns'
            coneIndices = opticDiskROI.indicesOfPointsInside(obj.coneRFpositionsMicrons);
        otherwise
            coneIndices = opticDiskROI.indicesOfPointsInside(obj.coneRFpositionsDegs);
    end
end
