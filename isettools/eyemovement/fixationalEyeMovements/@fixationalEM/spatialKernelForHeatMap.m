function [kernel, kernelSupport] = spatialKernelForHeatMap(obj)
% Create the spatial kernel for smoothing of the em heat map.
%
% Syntax:
%   [kernel, kernelSupport] = spatialKernalForHeatMap(obj)
%   [kernel, kernelSupport] = obj.spatialKernalForHeatMap()
%
% Description:
%    Create the spatial kernel with which to smooth the map of positions
%    visited during the recent em history. 
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
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
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