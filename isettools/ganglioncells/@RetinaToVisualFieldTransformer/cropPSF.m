function cropPSF(obj,maxSpatialSupportDegs)

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.theVlambdaWeightedPSFData.psfSupportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.theVlambdaWeightedPSFData.psfSupportYdegs) < maxSpatialSupportDegs);

    obj.theVlambdaWeightedPSFData.psfSupportXdegs = obj.theVlambdaWeightedPSFData.psfSupportXdegs(idx);
    obj.theVlambdaWeightedPSFData.psfSupportYdegs = obj.theVlambdaWeightedPSFData.psfSupportYdegs(idy);
    obj.theVlambdaWeightedPSFData.vLambdaWeightedData = obj.theVlambdaWeightedPSFData.vLambdaWeightedData(idy,idx);

    if (max(abs(maxSpatialSupportDegs(:))) > max(obj.theVlambdaWeightedPSFData.psfSupportXdegs))
        dx = obj.theVlambdaWeightedPSFData.psfSupportXdegs(2)-obj.theVlambdaWeightedPSFData.psfSupportXdegs(1);
        spatialSupportXdegs = dx:dx:max(abs(maxSpatialSupportDegs(:)));
        spatialSupportXdegs = [-fliplr(spatialSupportXdegs) 0 spatialSupportXdegs];
        obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs = spatialSupportXdegs;
        obj.theVlambdaWeightedPSFData.spatialSupportForRFmapYdegs = spatialSupportXdegs;
    else
        obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs = obj.theVlambdaWeightedPSFData.psfSupportXdegs;
        obj.theVlambdaWeightedPSFData.spatialSupportForRFmapYdegs = obj.theVlambdaWeightedPSFData.psfSupportXdegs;
    end

    fprintf('RetinaToVisualFieldTransformer: spatial support = [%2.2f ... +%2.2f] degs with %2.4f degs resolution.', ...
        min(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs), max(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs), dx);
       

    % Ensure we have unit volume
    obj.theVlambdaWeightedPSFData.vLambdaWeightedData = obj.theVlambdaWeightedPSFData.vLambdaWeightedData / sum(obj.theVlambdaWeightedPSFData.vLambdaWeightedData(:));
end