function setSizeToFOVForHexMosaic(obj,fov)
    obj.restoreOriginalResState();
    obj.setSizeToFOV(fov);
    obj.saveOriginalResState();
    obj.resampleGrid(obj.upSampleFactor);
end
