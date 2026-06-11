function currentPositionHeatLevel = updateHeatMap(obj, tStep)
% Update the heat map based on the eye movement position at tStep.
%
% Syntax:
%   currentPositionHeatLevel = updateHeatMap(obj, tStep)
%   currentPositionHeatLevel = obj.updateHeatMap(tStep)
%
% Description:
%    Update the heat map based on the eye movement position at the current
%    time step.
%
% Inputs:
%    obj                      - Object. A fixationalEM object.
%    tStep                    - Numeric. The current time step, in seconds.
%
% Outputs:
%    currentPositionHeatLevel - The heat level for the current position.
%
% Optional key/value pairs:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
%

% Re-center and scale to ArcMins
postStabilizationStep = 1 + tStep - obj.stabilizationStepsNum;
xPos = obj.scalarToArcMin * (obj.emPosTimeSeries(1, tStep) - ...
    obj.emPosTimeSeries(1, obj.stabilizationStepsNum));
yPos = obj.scalarToArcMin * (obj.emPosTimeSeries(2, tStep) - ...
    obj.emPosTimeSeries(2, obj.stabilizationStepsNum));

% Compute visited grid coords
gridNodesNum = length(obj.heatMapXYsupport);
node0 = (gridNodesNum - 1) / 2 + 1;
visitedGridCol = node0 + sign(xPos) * ...
    round(abs(xPos) / obj.heatMapSpatialSampleArcMin);
visitedGridRow = node0 + sign(yPos) * ...
    round(abs(yPos) / obj.heatMapSpatialSampleArcMin);

% At the beginning of the heatMapUpdateInterval, zero the spatial
% accumulation map
if (mod(postStabilizationStep - 1, obj.heatMapUpdateIntervalStepsNum) == 0)
    obj.emAccumMap = zeros(gridNodesNum, gridNodesNum);
end

% Accumulate the spatial map according to point visited
if ((visitedGridCol > 0) && (visitedGridCol <= gridNodesNum) && ...
        (visitedGridRow > 0) && (visitedGridRow <= gridNodesNum))
    obj.emAccumMap(visitedGridRow, visitedGridCol) = ...
        obj.emAccumMap(visitedGridRow, visitedGridCol) + 1;
end

% At the end of the heatMapUpdateInterval, perform spatial
% convolution of the accum map with the spatial kernel
if (mod(postStabilizationStep - 1, obj.heatMapUpdateIntervalStepsNum) ...
        == obj.heatMapUpdateIntervalStepsNum - 1)
    % Perform spatial convolution with heatMapSpatialKernel
    subSampledPostStabilizationStep = floor((postStabilizationStep - 1) ...
        / obj.heatMapUpdateIntervalStepsNum) + 1;
    obj.heatMapTimeSeriesIntermediate(subSampledPostStabilizationStep, ...
        :, :) = conv2(obj.emAccumMap, obj.heatMapSpatialKernel, 'same');

    % Perform temporal convolution with heatMapTemporalKernel
    tBinsForTemporalConvolution = 1:subSampledPostStabilizationStep;
    firstBinAffected = max([1 tBinsForTemporalConvolution(end) - ...
        size(obj.heatMapTemporalKernel, 1) + 1]);
    binsAffected = firstBinAffected : tBinsForTemporalConvolution(end);
    if (numel(binsAffected) == size(obj.heatMapTemporalKernel, 1))
        % Use the entire temporal filter
        instantaneousHeatMap = dot(obj.heatMapTimeSeriesIntermediate(...
            binsAffected, :, :), obj.heatMapTemporalKernel, 1);
    else
        % Use the valid portion of the temporal filter
        kernelPortion = ...
            obj.heatMapTemporalKernel(1:numel(binsAffected), :, :);
        instantaneousHeatMap = dot(obj.heatMapTimeSeriesIntermediate(...
            binsAffected, :, :), kernelPortion, 1);
    end
    % Inject copies of the instantaneousHeatMap map at appropriate delays
    % into the future.
    idx = postStabilizationStep + ...
        (0:obj.heatMapUpdateIntervalStepsNum - 1);
    obj.heatMapTimeSeries(idx, :, :) = ...
        repmat(instantaneousHeatMap, [numel(idx) 1 1]);
end

% Return heat level at current em position
currentPositionHeatLevel = 0;
if ((visitedGridCol>0) && (visitedGridCol <= gridNodesNum) && ...
    (visitedGridRow>0) && (visitedGridRow <= gridNodesNum))
        currentPositionHeatLevel = obj.heatMapTimeSeries(...
            postStabilizationStep, visitedGridRow, visitedGridCol);
end
% Keep a history
obj.currentPositionHeatLevelTimeSeries(1, postStabilizationStep) = ...
    currentPositionHeatLevel;

end