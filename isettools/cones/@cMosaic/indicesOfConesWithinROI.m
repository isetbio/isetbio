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

    % Validate the roi by calling the static method validateROI() of @cMosaic
    cMosaic.validateROI(roi);
    
    % Generate the roi outline by calling the static method generateOutline() of @cMosaic
    roiOutline = cMosaic.generateOutline(roi);
    
    % Convert roiOutline to microns
    if (strcmp(roi.units, 'degs'))
        roiOutlineMicrons = cMosaic.convertOutlineToMicrons(roiOutline,obj.micronsPerDegreeApproximation);
    else
        roiOutlineMicrons = roiOutline;
    end
    
    
    % Find indices of cones within the ROI border
    [in,on] = inpolygon( obj.coneRFpositionsMicrons(:,1), obj.coneRFpositionsMicrons(:,2),...
                         roiOutlineMicrons.x, roiOutlineMicrons.y);
    coneIndices = find((in|on));
    
end
