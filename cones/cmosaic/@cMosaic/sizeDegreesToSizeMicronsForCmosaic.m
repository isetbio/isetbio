% Method to transform a SIZE specified in units of visual degrees at a
% given eccentricity (also specified in degrees) to a SIZE in units of 
% retinal microns based on the @cMosaic configuration
function sizeMicrons = sizeDegreesToSizeMicronsForCmosaic(obj, sizeDegrees, eccentricityDegrees)

    if (~isempty(obj.micronsPerDegreeApproximation))
        % Linear transformation based on a constant
        sizeMicrons = sizeDegrees * obj.micronsPerDegreeApproximation;
    else
        if (~isempty(obj.customDegsToMMsConversionFunction))
            % Custom transformation based on arbitrary function
            if (numel(eccentricityDegrees) == 2) && (numel(sizeDegrees) > 1)
                % Treat eccentricityDegrees as a vector of (x,y) eccs
                eccentricityDegrees = sqrt(sum(eccentricityDegreess.^2,2));
            elseif (numel(eccentricityDegrees) == 1) || (numel(sizeDegrees) == 1)
                % Do nothing
            elseif (size(sizeDegrees,1) == size(eccentricityDegrees,1)) &&  (size(sizeDegrees,2)*size(eccentricityDegrees,2) == 1)
                % Do nothing
            elseif (size(sizeDegrees,2) == size(eccentricityDegrees,2)) &&  (size(sizeDegrees,1)*size(eccentricityDegrees,1) == 1)
                % Do nothing
            else
                error('dont know how to treat the passed dimensions');
            end
            
            leftCoordDegs  = (eccentricityDegrees - sizeDegrees/2);
            rightCoordDegs = (eccentricityDegrees + sizeDegrees/2);
            sizeMicrons =  1e3 * ( obj.customDegsToMMsConversionFunction(rightCoordDegs) ...
                                  -obj.customDegsToMMsConversionFunction(leftCoordDegs));               
        else
            % Watson's non-linear transformation
            sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegrees, eccentricityDegrees);
        end
    end
    
    
end
