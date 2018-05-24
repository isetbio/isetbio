function [kernel, kernelSupport] = spatialKernelForHeatMap(obj)
% Create the spatial kernel for the heat map.
%
% Syntax:
%   [kernel, kernelSupport] = spatialKernalForHeatMap(obj)
%   [kernel, kernelSupport] = obj.spatialKernalForHeatMap()
%
% Description:
%    Create the spatial kernel for the heat map.
%
% Inputs:
%    obj           - Object. A fixationalEM object.
%
% Outputs:
%    kernel        - Matrix. The spatial kernel heat map.
%    kernelSupport - Vector. The spatial kernal heat map's support.
%
% Optional key/value pairs:
%    None.
%

if (obj.heatMapKernelSpaceConstantArcMin < 0.01)
    sigma = 0.01;
else
    sigma = obj.heatMapKernelSpaceConstantArcMin;
end
halfPixelsNum = round(2.5 * sigma / obj.heatMapSpatialSampleArcMin);
ii = -halfPixelsNum:halfPixelsNum;
pixels = numel(ii);
kernel = fspecial('gaussian', pixels, sigma);
kernelSupport = ii * obj.heatMapSpatialSampleArcMin;

end