function cropPSF(obj,maxSpatialSupportDegs)

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.thePSFData.psfSupportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.thePSFData.psfSupportXdegs) < maxSpatialSupportDegs);

    obj.thePSFData.psfSupportXdegs = obj.thePSFData.psfSupportXdegs(idx);
    obj.thePSFData.psfSupportYdegs = obj.thePSFData.psfSupportXdegs(idy);
    obj.thePSFData.data = obj.thePSFData.data(idy,idx);

    if (max(abs(maxSpatialSupportDegs(:))) > max(obj.thePSFData.psfSupportXdegs))
        dx = obj.thePSFData.psfSupportXdegs(2)-obj.thePSFData.psfSupportXdegs(1);
        spatialSupportXdegs = dx:dx:max(abs(maxSpatialSupportDegs(:)));
        spatialSupportXdegs = [-fliplr(spatialSupportXdegs) 0 spatialSupportXdegs];
        obj.thePSFData.spatialSupportForRFmapXdegs = spatialSupportXdegs;
        obj.thePSFData.spatialSupportForRFmapYdegs = spatialSupportXdegs;
    else
        obj.thePSFData.spatialSupportForRFmapXdegs = obj.thePSFData.psfSupportXdegs;
        obj.thePSFData.spatialSupportForRFmapYdegs = obj.thePSFData.psfSupportXdegs;
    end
    
    % Ensure we have unit volume
    obj.thePSFData.data = obj.thePSFData.data / sum(obj.thePSFData.data(:));
end