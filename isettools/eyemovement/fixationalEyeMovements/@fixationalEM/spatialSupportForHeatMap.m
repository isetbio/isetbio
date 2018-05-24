function heatMapXYsupport = spatialSupportForHeatMap(obj)
% Method to compute the spatial support for the heatMap
%
% Syntax:
%   heatMapXYsupport = spatialSupportForHeatMap(obj)
%   heatMapXYsupport = obj.spatialSupportForHeatMap()
%
% Description:
%    Method to compute the spatial support for the heatMap
%
% Inputs:
%    obj              - Object. A fixationalEM object.
%
% Outputs:
%    heatMapXYsupport - Vector. The support vector for the fixational eye
%                       movement's heat map.
%
% Optional key/value pairs:
%    None.
%
gridNodesNum = ...
    ceil(obj.heatMapWidthArcMin / obj.heatMapSpatialSampleArcMin);
if (mod(gridNodesNum, 2) == 0), gridNodesNum = gridNodesNum + 1; end
node0 = (gridNodesNum-1)/2+1;
heatMapXYsupport = (-(node0 - 1):(node0 - 1)) * ...
    obj.heatMapSpatialSampleArcMin;

end