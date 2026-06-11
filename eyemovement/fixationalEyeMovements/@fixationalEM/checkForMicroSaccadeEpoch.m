function checkForMicroSaccadeEpoch(obj, tStep)
% Generate a microsaccade depending on the interval since the last one.
%
% Syntax:
%   obj.checkForMicroSaccadeEpoch(tStep)
%   checkForMicroSaccadeEpoch(obj, tStep)
%
% Description:
%    Check whether a microsaccade should be generated based on the inter-
%    saccade interval and statistics, and if so generate one according to
%    the strategy defined obj.microSaccadeType.
%
% Inputs:
%    obj   - Object. The fixationalEM object.
%    tStep - Numeric. The current time step
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments

% Compute probablility of generating a saccade base on interval since
% last saccade
intervalSinceLastMicroSaccadeSeconds = (tStep - ...
    obj.lastMicroSaccadeTimeStep) * obj.timeStepDurationSeconds;

if (round(intervalSinceLastMicroSaccadeSeconds * 1000) ~= ...
        obj.microSaccadeIntervalsMilliSecondsList(numel(...
        obj.microSaccadeOnsetStepIndices) + 1))
    % Not generating a microsaccade this time
    return;
end

% Ok, we are generating a saccade
%    Retrieve current position (aligned with respect to eye position at
%    the end of the stabilization time)
currentPosCentered = obj.emPosTimeSeries(:, tStep + 1) - ...
    obj.emPosTimeSeries(:, obj.stabilizationStepsNum);

% Generate a positional jitter vector
shape = 5;
scale = obj.microSaccadeTargetJitterArcMin / shape;
positionJitter = random('Gamma', shape, scale, [1 1]);
% ... randn(2, 1) .* obj.microSaccadeTargetJitterArcMin
positionJitterNorm = norm(positionJitter);

% Determine microsaccade target.
if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
    % Compute normalized currentHeatMap
    postStabilizationStep = 1 + tStep - obj.stabilizationStepsNum;
    currentHeatMap = squeeze(obj.heatMapTimeSeries(...
        postStabilizationStep, :, :));
    currentHeatMap = currentHeatMap / max(currentHeatMap(:));
    % Compute combined saccade target likelihood spatial map
    targetLikelihoodSpatialMap = obj.heatMapWeight * ...
        (1 - currentHeatMap) + (1 - obj.heatMapWeight) * ...
        obj.microSaccadeTargetLikelihoodSpatialMap;
    % Find point at which the targetLikelihoodSpatialMap is maximal
    targetPosArcMin = obj.posOfMaxLikelihood(...
        targetLikelihoodSpatialMap);
    targetPos = targetPosArcMin / obj.scalarToArcMin;
    % Compute the displacement vector
    microSaccadeJumpArcMin = (targetPos' - currentPosCentered) * ...
        obj.scalarToArcMin + positionJitter;
    % Compute the duration based on the displacement vector and a speed
    % from the obj.microSaccadeSpeedDegsPerSecondList
    saccadeIndex = numel(obj.microSaccadeOnsetStepIndices) + 1;
    speedArcMinPerMillisecond = ...
        obj.microSaccadeSpeedDegsPerSecondList(saccadeIndex) / 1000 * 60;
    durationMilliSeconds = round(norm(microSaccadeJumpArcMin) / ...
        speedArcMinPerMillisecond);
    durationMilliSeconds = max([...
        obj.microSaccadeMinDurationMilliSeconds durationMilliSeconds]);

elseif (strcmp(obj.microSaccadeType, 'stats based'))
    % Decide whether we will do a corrective saccade based on how far
    % away we are
    correctionDecisionFactor = 1 - exp(-(norm(currentPosCentered * ...
        obj.scalarToArcMin) / ...
        (0.5 * obj.fixationMapSpaceConstantArcMin)) ^ 2);
    doCorrectiveSaccade = correctionDecisionFactor > rand(1, 1);
    if (doCorrectiveSaccade)
        % Corrective saccade
        % direction toward the fixation point
        theta = atan2(currentPosCentered(2), currentPosCentered(1)) + ...
            pi + randn(1, 1) * obj.microSaccadeDirectionJitterDegs / ...
            180 * pi;

        % Distance to fixation
        distanceToFixationPointArcMin = ...
            norm(squeeze(currentPosCentered)) * obj.scalarToArcMin;
        % Add some jitter
        distanceToFixationPointArcMin = ...
            distanceToFixationPointArcMin + positionJitterNorm;

        % Find the index of closest amplitude from the
        % obj.microSaccadeAmplitudesArcMinList
        [~, idx] = min(abs(distanceToFixationPointArcMin - ...
            obj.microSaccadeAmplitudesArcMinList(...
            obj.validMicrosaccadeAmplitudeIndices)));
        idx = obj.validMicrosaccadeAmplitudeIndices(idx);
    else
        % Random saccade
        %    random direction, random amplitude, but one that does not
        %    take us further away than some threshold distance
        theta = rand(1, 1) * 2 * pi;

        % Find an amplitude index that is within the
        % obj.fixationMapSpaceConstantArcMin region
        possibleJumps = (obj.microSaccadeAmplitudesArcMinList(...
            obj.validMicrosaccadeAmplitudeIndices))' * ...
            [cos(theta) sin(theta)];
        possibleEndPoints = bsxfun(@plus, possibleJumps, ...
            (currentPosCentered') * obj.scalarToArcMin);
        possibleDistancesToFixationPoint = sqrt(sum(...
            possibleEndPoints .^ 2, 2));
        idx = find(possibleDistancesToFixationPoint < ...
            obj.fixationMapSpaceConstantArcMin);
        if (isempty(idx))
            [~, idx] = min(possibleDistancesToFixationPoint);
        else
            randomIdx = randperm(numel(idx));
            idx = idx(randomIdx(1));
        end
        idx = obj.validMicrosaccadeAmplitudeIndices(idx);
    end

    % Retrieve the amplitude at the computed index
    amplitudeArcMin = obj.microSaccadeAmplitudesArcMinList(idx);

    % Retrieve the correspondng duration
    durationMilliSeconds = ...
        obj.microSaccadeDurationsMillisecondsList(idx);

    % Remove this amplitude from the list of valid indices (sampling
    % without replacement)
    obj.validMicrosaccadeAmplitudeIndices = setdiff(...
        obj.validMicrosaccadeAmplitudeIndices, idx);

    % Compute the displacement vector
    microSaccadeJumpArcMin = amplitudeArcMin * [cos(theta); sin(theta)];
else
    error('Unknown microsaccade type: ''%s''\n', obj.microSaccadeType);
end

% Saccade jump in native units
saccadeJump = microSaccadeJumpArcMin / obj.scalarToArcMin;
obj.microSaccadeResidualPath = ones(2, durationMilliSeconds);
obj.microSaccadeResidualPath(1, :) = saccadeJump(1, 1) / ...
    durationMilliSeconds;
obj.microSaccadeResidualPath(2, :) = saccadeJump(2, 1) / ...
    durationMilliSeconds;

% Insert the saccade jump
endingPosition = obj.emPosTimeSeries(:, tStep+1) + saccadeJump;
obj.emPosTimeSeries(:, tStep + 1) = ...
    obj.emPosTimeSeries(:, tStep + 1) + obj.microSaccadeResidualPath(:, 1);
% Update the microSaccadeResidualPath
obj.updateMicroSaccadeResidualPath();

% Keep a record of the saccade onset time indices, saccade jumps, and
% saccade intervals
obj.microSaccadeOnsetStepIndices(...
    numel(obj.microSaccadeOnsetStepIndices) + 1) = tStep + 1;
saccadeIndex = numel(obj.microSaccadeOnsetStepIndices);
obj.microSaccadeAmplitudesArcMin(saccadeIndex) = ...
    norm(microSaccadeJumpArcMin);
obj.microSaccadeSpeedsDegsPerSecond(saccadeIndex) = ...
    (norm(microSaccadeJumpArcMin) / 60) / (durationMilliSeconds / 1000);
obj.microSaccadeIntervalsMilliseconds(saccadeIndex) = ...
    intervalSinceLastMicroSaccadeSeconds * 1000;
if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
    obj.microSaccadeTargetLikelihoodMapsList(saccadeIndex, :, :) = ...
        targetLikelihoodSpatialMap;
end
obj.microSaccadePositions(saccadeIndex, :) = obj.scalarToArcMin * ...
    (endingPosition - obj.emPosTimeSeries(:, obj.stabilizationStepsNum))';

% Update time of last saccade
obj.lastMicroSaccadeTimeStep = obj.microSaccadeOnsetStepIndices(end);

end