function removeConesWithinOpticNerveHead(obj)
% Remove cones located within the optic disk
%
% Syntax:
%   obj.removeConesWithinOpticNerveHead()
%
% Description:
%    Remove cones located within the optic disk
%
% Inputs:
%    obj                 - A @cMosaic object
%
% Outputs:                 None

    % Convert OD center from degs to microns
    if (~isempty(obj.micronsPerDegreeApproximation))
        % Convert outline from degs to microns using the passed
        % microns/deg approximation
        centerMicrons = obj.opticDisk.centerDegs * obj.micronsPerDegreeApproximation;
    else
        centerMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(obj.opticDisk.centerDegs);
    end
        
    % The optic disk is located in the nasal retinas. In terms of the RE
    % visual field, this corresponds to a negative x-coord for the left eye
    % and a positive x-coord for the right eye
    if (strcmp(obj.whichEye, 'left eye'))
        centerMicrons(1) = -centerMicrons(1);
    end
    % y-position is inversed in the visual field of both eye
    centerMicrons(2) = centerMicrons(2);
    
    % Construct optic disk ROI struct in microns
    opticDiskROI = struct(...
        'units', 'microns', ...
        'shape', 'ellipse', ...
        'center', centerMicrons, ...
        'minorAxisDiameter', obj.opticDisk.horizontalDiameterMM*1e3, ...
        'majorAxisDiameter', obj.opticDisk.verticalDiameterMM*1e3, ...
        'rotation', 0);
    
    % Find indices of cones lying inside the optic disk ellipsoid
    idxInside = obj.indicesOfConesWithinROI(opticDiskROI);
    
    % Find indices of cones outside the optic disk
    idxOutside = setdiff(1:size(obj.coneRFpositionsDegs,1), idxInside);

    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idxOutside,:);
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idxOutside,:);
    obj.coneRFspacingsDegs = obj.coneRFspacingsDegs(idxOutside);
    obj.coneRFspacingsMicrons = obj.coneRFspacingsMicrons(idxOutside);
end
                  