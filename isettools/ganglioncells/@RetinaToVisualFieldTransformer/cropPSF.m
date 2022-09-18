function cropPSF(obj,maxSpatialSupportDegs)

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.thePSFData.supportXdegs) < maxSpatialSupportDegs);
    idy = find(abs(obj.thePSFData.supportYdegs) < maxSpatialSupportDegs);

    if (maxSpatialSupportDegs > max(obj.thePSFData.supportXdegs(:)))
        fprintf(2, '>>>>cropPSF:: PSF must be zero padded to cover the needed spatial support');
        pause;
    end

    obj.thePSFData.supportXdegs = obj.thePSFData.supportXdegs(idx);
    obj.thePSFData.supportYdegs = obj.thePSFData.supportYdegs(idy);
    obj.thePSFData.data = obj.thePSFData.data(idy,idx);
    
    % Ensure we have unit volume
    obj.thePSFData.data = obj.thePSFData.data / sum(obj.thePSFData.data(:));
end