function trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds)
% Perform a trim to only keep samples after the stabilization time
%
% Syntax:
%   trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds)
%   obj.trimRecenterAndResampleTimeSeries(sampleDurationSeconds)
%
% Description:
%    Trim the data and resample in such a fashion that only the samples
%    after the stabilization time are kept.
%
% Inputs:
%    obj                   - Object. A fixationalEM object.
%    sampleDurationSeconds - Numeric. The sample duration, in seconds.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Trim: only keep samples after the stabilization time
keptSteps = obj.stabilizationStepsNum + 1:obj.tStepsNum - 1;
obj.timeAxis = obj.timeAxis(keptSteps);
obj.emPosTimeSeries = obj.emPosTimeSeries(:, keptSteps);
obj.heatMapTimeSeries = obj.heatMapTimeSeries(1:numel(keptSteps), :, :);
obj.velocityTimeSeries = obj.velocityTimeSeries(keptSteps);
obj.microSaccadeOnsetStepIndices = obj.microSaccadeOnsetStepIndices - ...
    obj.stabilizationStepsNum;
obj.microSaccadeOnsetStepIndices = obj.microSaccadeOnsetStepIndices(...
    obj.microSaccadeOnsetStepIndices >= 0);

% Re-center
obj.timeAxis = obj.timeAxis - obj.timeAxis(1);
obj.emPosTimeSeries(1, :) = ...
    obj.emPosTimeSeries(1, :) - obj.emPosTimeSeries(1, 1);
obj.emPosTimeSeries(2, :) = ...
    obj.emPosTimeSeries(2, :) - obj.emPosTimeSeries(2, 1);

% Resample in time according to passed sampleDurationSeconds
if (abs(obj.timeStepDurationSeconds - sampleDurationSeconds) > ...
        100 * eps(sampleDurationSeconds))
    % Resampled time axis
    oldTimeAxis = obj.timeAxis;
    obj.timeAxis = 0:sampleDurationSeconds:obj.timeAxis(end);
    obj.emPosTimeSeries = obj.smartInterpolation(...
        oldTimeAxis, obj.emPosTimeSeries, obj.timeAxis);
    obj.velocityTimeSeries = obj.smartInterpolation(...
        oldTimeAxis, obj.velocityTimeSeries, obj.timeAxis);
    obj.heatMapTimeSeries = obj.smartInterpolation(...
        oldTimeAxis, obj.heatMapTimeSeries, obj.timeAxis);
end

end
