% Method to compute the spatial support for the heatMap
function heatMapXYsupport = spatialSupportForHeatMap(obj)
    gridNodesNum = ceil(obj.heatMapWidthArcMin / obj.heatMapSpatialSampleArcMin);
    if (mod(gridNodesNum,2) == 0)
        gridNodesNum = gridNodesNum + 1;
    end
    node0 = (gridNodesNum-1)/2+1;
    heatMapXYsupport = (-(node0-1):(node0-1))*obj.heatMapSpatialSampleArcMin;
end