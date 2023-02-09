function cropPSF(obj, maxSpatialSupportDegs)


    % Find indices of samples within the maxSpatialSupport
    idx = find(abs(obj.spectrallyWeightedPSFData.psfSupportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.spectrallyWeightedPSFData.psfSupportYdegs) < maxSpatialSupportDegs);

    % Crop the spatial support
    obj.spectrallyWeightedPSFData.psfSupportXdegs = ...
        obj.spectrallyWeightedPSFData.psfSupportXdegs(idx);

    obj.spectrallyWeightedPSFData.psfSupportYdegs = ...
        obj.spectrallyWeightedPSFData.psfSupportYdegs(idy);

    obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs = ...
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs(idx);

    obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs = ...
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs (idy);


    % Crop all weighted PSFs
    obj.spectrallyWeightedPSFData.LconeWeighted = ...
        obj.spectrallyWeightedPSFData.LconeWeighted(idy,idx);

    obj.spectrallyWeightedPSFData.MconeWeighted = ...
        obj.spectrallyWeightedPSFData.MconeWeighted(idy,idx);

    obj.spectrallyWeightedPSFData.LMconeWeighted = ...
        obj.spectrallyWeightedPSFData.LMconeWeighted(idy,idx);

    % Ensure we have unit volume
    obj.spectrallyWeightedPSFData.LconeWeighted  = ...
        obj.spectrallyWeightedPSFData.LconeWeighted / sum(obj.spectrallyWeightedPSFData.LconeWeighted(:));

    obj.spectrallyWeightedPSFData.MconeWeighted  = ...
        obj.spectrallyWeightedPSFData.MconeWeighted / sum(obj.spectrallyWeightedPSFData.MconeWeighted(:));

    obj.spectrallyWeightedPSFData.LMconeWeighted  = ...
        obj.spectrallyWeightedPSFData.LMconeWeighted / sum(obj.spectrallyWeightedPSFData.LMconeWeighted(:));

    % Increase spatial support for RFmaps if needed
    if (max(abs(maxSpatialSupportDegs(:))) > 1.02 * max(abs(obj.spectrallyWeightedPSFData.psfSupportXdegs(:))))
        
        dx = obj.spectrallyWeightedPSFData.psfSupportXdegs(2)-obj.spectrallyWeightedPSFData.psfSupportXdegs(1);
        spatialSupportXdegs = dx:dx:max(abs(maxSpatialSupportDegs(:)));
        spatialSupportXdegs = [-fliplr(spatialSupportXdegs) 0 spatialSupportXdegs];
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs = spatialSupportXdegs;
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs = spatialSupportXdegs;

        fprintf('Sspatial support for PSF [%2.4f .. %2.4f] degs\n', ...
            min(obj.spectrallyWeightedPSFData.psfSupportXdegs), max(obj.spectrallyWeightedPSFData.psfSupportXdegs));
        fprintf('Increasing spatial support for RF maps: [%2.4f .. %2.4f] degs\n', ...
            min(obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs ), max(obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs ));

    end

end
