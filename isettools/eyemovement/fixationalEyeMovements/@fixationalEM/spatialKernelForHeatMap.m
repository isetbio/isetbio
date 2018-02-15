function [kernel, kernelSupport] = spatialKernelForHeatMap(obj)
    if (obj.heatMapKernelSpaceConstantArcMin < 0.01)
        sigma = 0.01;
    else
        sigma = obj.heatMapKernelSpaceConstantArcMin;
    end
    halfPixelsNum = round(2.5*sigma/obj.heatMapSpatialSampleArcMin);
    ii = -halfPixelsNum:halfPixelsNum;
    pixels = numel(ii);
    kernel = fspecial('gaussian', pixels, sigma);
    kernelSupport = ii*obj.heatMapSpatialSampleArcMin;
end