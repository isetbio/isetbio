function [odStructMicrons, odStructDegs] = odStruct(obj)
    % Convert OD center from degs to microns
    centerMicrons = obj.distanceDegreesToDistanceMicronsForCmosaic(obj.opticDisk.centerDegs);
        
    % The optic disk is located in the nasal retinas. In terms of the RE
    % visual field, this corresponds to a negative x-coord for the left eye
    % and a positive x-coord for the right eye
    if (strcmp(obj.whichEye, 'left eye'))
        centerMicrons(1) = -centerMicrons(1);
    end
    % y-position is inversed in the visual field of both eye
    centerMicrons(2) = -centerMicrons(2);

    % Construct optic disk ROI struct in microns
    odStructMicrons = struct(...
        'units', 'microns', ...
        'shape', 'ellipse', ...
        'center', centerMicrons, ...
        'minorAxisDiameter', obj.opticDisk.horizontalDiameterMM*1e3, ...
        'majorAxisDiameter', obj.opticDisk.verticalDiameterMM*1e3, ...
        'rotation', obj.opticDisk.rotationDegs);
    
    % Construct optic disk ROI struct in degs
    % Convert OD size from mm to degs
    minorAxisDiameterDegs = obj.distanceMicronsToDistanceDegreesForCmosaic(1e3 * obj.opticDisk.horizontalDiameterMM);
    majorAxisDiameterDegs = obj.distanceMicronsToDistanceDegreesForCmosaic(1e3 * obj.opticDisk.verticalDiameterMM);
    
    centerDegs = obj.opticDisk.centerDegs(1);
    if (strcmp(obj.whichEye, 'left eye'))
        centerDegs(1) = -obj.opticDisk.centerDegs(1);
    end
    % y-position is inversed in the visual field of both eye
    centerDegs(2) = -obj.opticDisk.centerDegs(2);
    
    odStructDegs = struct(...
        'units', 'degs', ...
        'shape', 'ellipse', ...
        'center', centerDegs, ...
        'minorAxisDiameter', minorAxisDiameterDegs, ...
        'majorAxisDiameter', majorAxisDiameterDegs, ...
        'rotation',obj.opticDisk.rotationDegs);
end

    