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

    % Validate the roi
    validateROI(roi);
    
    % Generate the roi outline
    roiOutline = generateOutline(roi);
    
    % Convert roiOutline to microns
    if (strcmp(roi.units, 'degs'))
        roiOutlineMicrons = convertOutlineToMicrons(roiOutline,obj);
    else
        roiOutlineMicrons = roiOutline;
    end
    
    
    % Find indices of cones within the ROI border
    [in,on] = inpolygon( obj.coneRFpositionsMicrons(:,1), obj.coneRFpositionsMicrons(:,2),...
                         roiOutlineMicrons.x, roiOutlineMicrons.y);
    coneIndices = find((in|on));
    
end

function validateROI(roi)
    % roi must be a struct
    assert(isstruct(roi), 'Region of interest must be a struct.');
    
    % roi must have a 'shape' field which is set to either 'rect' or 'ellipse'
    assert(isfield(roi, 'shape'), 'roi struct must have a ''shape'' field.');
    assert(ismember(roi.shape, {'rect', 'ellipse'}),'roi.shape can be either ''rect'' or ''ellipse''.');
    
    % roi must have a 'units' field which is set to either 'degs' or 'microns'
    assert(isfield(roi, 'units'), 'roi struct must have a ''units'' field.');
    assert(ismember(roi.units, {'degs', 'microns'}),'roi.units can be either ''degs'' or ''microns''.');
    
    % roi must have a 'center' field
    assert(isfield(roi, 'center'), 'roi struct must have a ''center'' field.');
    
    if (strcmp(roi.shape, 'rect'))
        % rect roi must have 'width' and 'height' fields
        assert(isfield(roi, 'width'), 'roi struct must have a ''width'' field.');
        assert(isfield(roi, 'height'), 'roi struct must have a ''height'' field.');
    else
        % ellipse roi must have 'rotation', 'minorAxisDiameter' and 'majorAxisDiameter' fields
        assert(isfield(roi, 'rotation'), 'roi struct must have a ''rotation'' field.');
        assert(isfield(roi, 'minorAxisDiameter'), 'roi struct must have a ''minorAxisDiameter'' field.');
        assert(isfield(roi, 'minorAxisDiameter'), 'roi struct must have a ''majorAxisDiameter'' field.');
    end
end

function roiOutline = generateOutline(roi)
    if (strcmp(roi.shape, 'rect'))
        roiOutline.x = roi.center(1) + 0.5*roi.width  * [-1 -1  1  1 -1];
        roiOutline.y = roi.center(2) + 0.5*roi.height * [-1  1  1 -1 -1];
    else
        a = 0.5*roi.minorAxisDiameter;
        b = 0.5*roi.majorAxisDiameter;
        cosOutline = cosd(0:10:360);
        sinOutline = sind(0:10:360);
        roiOutline.x = roi.center(1) + a*cosOutline*cosd(roi.rotation) - b*sinOutline*sind(roi.rotation);
        roiOutline.y = roi.center(2) + a*cosOutline*sind(roi.rotation) + b*sinOutline*cosd(roi.rotation);
    end
end

function roiOutlineMicrons = convertOutlineToMicrons(roiOutlineDegs,obj)
    % Convert degs to microns
    if (~isempty(obj.micronsPerDegreeApproximation))
        % Convert outline from degs to microns using the passed microns/deg approximation
        roiOutlineMicrons.x = roiOutlineDegs.x * obj.micronsPerDegreeApproximation;
        roiOutlineMicrons.y = roiOutlineDegs.y * obj.micronsPerDegreeApproximation;
    else
        roiOutlineMicrons.x = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(roiOutlineDegs.x); 
        roiOutlineMicrons.y = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(roiOutlineDegs.y); 
    end
end
