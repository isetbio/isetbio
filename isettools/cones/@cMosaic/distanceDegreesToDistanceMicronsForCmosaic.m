% Method to transform a DISTANCE specified in units of retinal
% microns to DISTANCE in units of visual degrees based on the @cMosaic configuration
function microns = distanceDegreesToDistanceMicronsForCmosaic(obj, degrees)

    if (~isempty(obj.micronsPerDegreeApproximation))
        % Linear transformation based on a constant
        microns = degrees * obj.micronsPerDegreeApproximation;
    else
        if (~isempty(obj.customDegsToMMsConversionFunction))
            % Custom transformation based on arbitrary function
            microns = 1e3 * obj.customDegsToMMsConversionFunction(degrees);
        else
            % Watson's non-linear transformation
            microns = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(degrees);
        end
    end

end

