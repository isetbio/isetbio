function cropPSF(obj,maxSpatialSupportDegs)

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.thePSFData.supportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.thePSFData.supportYdegs) < maxSpatialSupportDegs);

    obj.thePSFData.supportXdegs = obj.thePSFData.supportXdegs(idx);
    obj.thePSFData.supportYdegs = obj.thePSFData.supportYdegs(idy);
    obj.thePSFData.data = obj.thePSFData.data(idy,idx);
    
    obj.theCircularPSFData.supportXdegs = obj.theCircularPSFData.supportXdegs(idx);
    obj.theCircularPSFData.supportYdegs = obj.theCircularPSFData.supportYdegs(idy);
    obj.theCircularPSFData.data = obj.theCircularPSFData.data(idy,idx);
end