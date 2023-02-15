function cropPSF(obj,maxSpatialSupportDegs)

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.theSpectrallyWeightedPSFData.psfSupportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.theSpectrallyWeightedPSFData.psfSupportYdegs) < maxSpatialSupportDegs);

    % Crop the spatial support
    obj.theSpectrallyWeightedPSFData.psfSupportXdegs = obj.theSpectrallyWeightedPSFData.psfSupportXdegs(idx);
    obj.theSpectrallyWeightedPSFData.psfSupportYdegs = obj.theSpectrallyWeightedPSFData.psfSupportYdegs(idy);

    % Crop all weighted PSFs
    obj.theSpectrallyWeightedPSFData.LconeWeighted = ...
        obj.theSpectrallyWeightedPSFData.LconeWeighted(idy,idx);

    obj.theSpectrallyWeightedPSFData.MconeWeighted = ...
        obj.theSpectrallyWeightedPSFData.MconeWeighted(idy,idx);

    obj.theSpectrallyWeightedPSFData.LMconeWeighted = ...
        obj.theSpectrallyWeightedPSFData.LMconeWeighted(idy,idx);


    if (max(abs(maxSpatialSupportDegs(:))) > max(obj.theSpectrallyWeightedPSFData.psfSupportXdegs))
        dx = obj.theSpectrallyWeightedPSFData.psfSupportXdegs(2)-obj.theSpectrallyWeightedPSFData.psfSupportXdegs(1);
        spatialSupportXdegs = dx:dx:max(abs(maxSpatialSupportDegs(:)));
        spatialSupportXdegs = [-fliplr(spatialSupportXdegs) 0 spatialSupportXdegs];
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs = spatialSupportXdegs;
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapYdegs = spatialSupportXdegs;
    else
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs = obj.theSpectrallyWeightedPSFData.psfSupportXdegs;
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapYdegs = obj.theSpectrallyWeightedPSFData.psfSupportXdegs;
    end

    fprintf('RetinaToVisualFieldTransformer: spatial support = [%2.2f ... +%2.2f] degs with %2.4f degs resolution.', ...
        min(obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs), max(obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs), dx);
       
    % Ensure we have unit volume
    obj.theSpectrallyWeightedPSFData.LconeWeighted  = ...
        obj.theSpectrallyWeightedPSFData.LconeWeighted / sum(obj.theSpectrallyWeightedPSFData.LconeWeighted(:));

    obj.theSpectrallyWeightedPSFData.MconeWeighted  = ...
        obj.theSpectrallyWeightedPSFData.MconeWeighted / sum(obj.theSpectrallyWeightedPSFData.MconeWeighted(:));

    obj.theSpectrallyWeightedPSFData.LMconeWeighted  = ...
        obj.theSpectrallyWeightedPSFData.LMconeWeighted / sum(obj.theSpectrallyWeightedPSFData.LMconeWeighted(:));
end