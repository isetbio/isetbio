% Method to transform a DISTANCE specified in units of visual
% degrees to DISTANCE in units of retinal microns based on the @cMosaic configuration
function degrees = distanceMicronsToDistanceDegreesForCmosaic(obj, microns)

    if (~isempty(obj.micronsPerDegreeApproximation))
        % Linear transformation based on a constant
        degrees = microns / obj.micronsPerDegreeApproximation;
    else
        if (~isempty(obj.customMMsToDegsConversionFunction))
            % Custom transformation based on arbitrary function
            degrees = obj.customMMsToDegsConversionFunction(microns*1e-3);
        else
            % Watson's non-linear transformation
            degrees = RGCmodels.Watson.convert.rhoMMsToDegs(microns*1e-3);
        end
    end
end


     
        