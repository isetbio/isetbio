% Method to transform a SIZE specified in units of retinal microns at a
% given eccentricity (also specified in microns) to a SIZE in units of 
% visual degrees based on the @cMosaic configuration
function sizeDegrees = sizeMicronsToSizeDegreesForCmosaic(obj, sizeMicrons, eccentricityMicrons)
    if (~isempty(obj.micronsPerDegreeApproximation))
        % Linear transformation based on a constant
        sizeDegrees = sizeMicrons / obj.micronsPerDegreeApproximation;
    else
        if (~isempty(obj.customMMsToDegsConversionFunction))
            % Custom transformation based on arbitrary function
            if (numel(eccentricityMicrons) == 2) && (numel(sizeMicrons) > 1)
                % Treat eccentricityMicrons as a vector of (x,y) eccs
                eccentricityMicrons = sqrt(sum(eccentricityMicrons.^2,2));
            elseif (numel(eccentricityMicrons) == 1) || (numel(sizeMicrons) == 1)
                % Do nothing
            elseif (size(sizeMicrons,1) == size(eccentricityMicrons,1)) &&  (size(sizeMicrons,2)*size(eccentricityMicrons,2) == 1)
                % Do nothing
            elseif (size(sizeMicrons,2) == size(eccentricityMicrons,2)) &&  (size(sizeMicrons,1)*size(eccentricityMicrons,1) == 1)
                % Do nothing
            else
                error('dont know how to treat the passed dimensions');
            end
            
            leftCoordMM  = (eccentricityMicrons - sizeMicrons/2)*1e-3;
            rightCoordMM = (eccentricityMicrons + sizeMicrons/2)*1e-3;
            sizeDegrees = obj.customMMsToDegsConversionFunction(rightCoordMM) - ...
                          obj.customMMsToDegsConversionFunction(leftCoordMM);
        else
            % Watson's non-linear transformation
            sizeDegrees = RGCmodels.Watson.convert.sizeRetinalMicronsToSizeVisualDegs(sizeMicrons, eccentricityMicrons);
        end
    end
    
end
