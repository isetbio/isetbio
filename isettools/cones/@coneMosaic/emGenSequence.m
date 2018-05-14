function fixEMobj = emGenSequence(obj, nEyeMovements, varargin)
% Generate sequence of fixational eye movements for a rect cone mosaic
%
% Syntax:
%	[nTrialsPos, nTrialsPosMicrons] = emGenSequence(...
%       obj, nEyeMovements, [varargin])
%
% Description:
%    The eye movement samples are created at the same temporal sample rate
%    as the cone integration time. We only update the position at the
%    beginning of each integration time.
%
%
% Inputs:
%     obj               - rect cone mosaic object
%     nEyeMovements     - number of eye movements to generate
%
% Ouputs:
%     nTrialsPos        - eye positions in a tensor of size 
%                            nTrials x nEyeMovements x 2 matrix.  
%                         Units are cone positions
%     nTrialsPosMicrons - eye positions in a tensor of size
%                            nTrials x nEyeMovements x 2 matrix
%                         Units are microns.
%
% Optional key/value pairs:
%    'em'               - Fixational eye movement structure, see
%                           fixationalEM for details
%    'rSeed'            - Random seed to be used
%    'nTrials'          - Multiple trial case, default = 1
%
% See Also:
%     fixationalEM
%

% History:
%    xx/xx/16  HJ/BW    ISETBIO Team, 2016
%    11/06/17  ncp      Added line to make drift magnitude independent of
%                       sample time.
%    11/06/17  dhb/npc  Added comments on microsaccade algorithm.
%    11/07/17  dhb      More cleaning and robustness.
%    02/26/18  jnm      Formatting, fix example
%    05/13/18  baw      Re-wrote for fixationalEM class.

% Examples:
%{
 scene = sceneCreate('mackay'); scene = sceneSet(scene,'fov',1);
 oi = oiCreate; oi = oiCompute(oi,scene);
 cm = coneMosaic; 

 nEyeMoves = 50;
 cm.emGenSequence(nEyeMoves, 'nTrials', 1);
 cm.compute(oi);  cm.window;

 cm.emGenSequence(nEyeMoves, 'nTrials', 1, 'microsaccade type','heatmap/fixation based');
 cm.compute(oi);  cm.window;

%}

%% parse inputs
p = inputParser;
varargin = ieParamFormat(varargin);

p.addRequired('obj', @(x)(isa(x,'coneMosaic')));
p.addRequired('nEyeMovements', @isscalar);

% Either create the default, or the user creates it with special parameters
% and provides it
p.addParameter('em', fixationalEM, @(x)(isa(x,'fixationalEM')));

% User can set the microSaccadeType, but nothing else, about the
% fixationalEM.
validTypes = {'none','stats based','heatmap/fixation based'};
p.addParameter('microsaccadetype', 'none', @(x)(ismember(x,validTypes)));
p.addParameter('ntrials', 1, @isscalar);
p.addParameter('rseed', 1, @isscalar);
p.addParameter('computevelocity', [], @islogical);
p.addParameter('useparfor', false, @islogical);

% set parameters
p.parse(obj, nEyeMovements, varargin{:});

fixEMobj   = p.Results.em;
nTrials    = p.Results.ntrials;
randomSeed = p.Results.rseed;
if ~isempty(p.Results.rseed), rng(p.Results.rseed); end

microSaccadeType = p.Results.microsaccadetype;

%% Start the calculation

fixEMobj.microSaccadeType = microSaccadeType;
fixEMobj.randomSeed = randomSeed;
fixEMobj.computeForConeMosaic(obj,nEyeMovements, ...
    'nTrials',nTrials, ...
    'rSeed', randomSeed);

%% Set the cone eye movement positions variable

obj.emPositions = fixEMobj.emPos;

end
%{
% Deprecated calculation from HJ days.

for nn=1:nTrials
    %% generate eye movement for tremor
    if emFlag(1) || microSaccadesOnlyFlag 
        % Load parameters
        amplitude = emGet(em, 'tremor amplitude', 'cones/sample', params);
        interval = emGet(em, 'tremor interval');
        intervalSD = emGet(em, 'tremor interval SD');

        % Compute time of tremor occurs
        t = interval + randn(nEyeMovements, 1) * intervalSD;
        t(t < 0.001) = 0.001; % get rid of negative values
        tPos = cumsum(t);
        tPos = max(1, round(tPos / sampTime));
        indx = 1:find(tPos <= nEyeMovements, 1, 'last');
        tPos = tPos(indx);

        % Generate random step at the selected times
        direction = rand(length(tPos), 1);

        % Unit length direction
        pos(tPos, :) = amplitude * [direction sqrt(1 - direction .^ 2)];

        pos(tPos, :) = bsxfun(@times, pos(tPos, :), t(indx(:)) / sampTime);
        pos = pos .* (2 * (randn(size(pos)) > 0) - 1); % shuffle the sign
        pos = cumsum(pos, 1);
    end

    % Generate eye movement for drift
    if emFlag(2) || microSaccadesOnlyFlag  
        % Load Parameters
        speed = emGet(em, 'drift speed', 'cones/sample', params);
        speedSD = emGet(em, 'drift speed SD', 'cones/sample', params);

        % Make drift magnitude independent of sample time.
        correctionForSampleTime = 1000 * sampTime;

        % Generate random move at each sample time
        theta = 360 * randn + 0.1 * correctionForSampleTime * ...
            (1 : nEyeMovements)';
        direction = [cosd(theta) sind(theta)];
        s = speed + speedSD * randn(nEyeMovements, 1);
        pos = filter(1, [1 -1], bsxfun(@times, direction, s)) + pos;
    end

    if (microSaccadesOnlyFlag)
        % Save eye movements due to drift and tremor
        driftTremorEMpos = pos;
    end

    % generate eye movement for micro-saccade
    if emFlag(3)
        % Load microsaccade model parameters
        interval = emGet(em, 'msaccade interval');
        intervalSD = emGet(em, 'msaccade interval SD');
        dirSD = emGet(em, 'msaccade dir SD', 'deg');
        speed = emGet(em, 'msaccade speed', 'cones/sample', params);
        speedSD = emGet(em, 'msaccade speed SD', 'cones/sample', params);

        % For the logic of this routine to be correct, we are assuming
        % that the integration time (sampTime) is one second or less.
        % Check that condition here.
        if (sampTime > 1)
            error(['Sample time (aka integration time is too long ' ...
                '(greater than 1 sec)']);
        end

        % Compute microsaccade occurence times.
        %
        % The use of randn means that we might end up with some
        % negative intervals, we just get rid of these by substituting
        % in a small interval
        minInterval = 0.010;
        saccadeIntervals = interval + randn(nEyeMovements, 1) * intervalSD;
        saccadeIntervals(saccadeIntervals < 0) = minInterval;

        % OLD: Old way of getting rid of negative intervals, seems like it
        % will produce a buildup of intervals around a hard coded 0.3, and
        % that it could still produce negative intervals.
        % t(t < 0.3) = 0.3 + 0.1 * rand;

        % Create list of discrete sample times
        discreteSampleTimes = (0:(nEyeMovements - 1)) * sampTime;

        % Ensure at least one saccade
        if (isempty(saccadeIntervals))
            saccadeIntervals = 0.5 * nEyeMovements / sampTime;
        end

        % Convert saccade intervals to saccade times, and don't go past the
        % end of the sample times. We need to do the latter check because
        % the list of ossible saccade times was created to be very long, at
        % least for reasonable integration times.
        saccadeTimes = cumsum(saccadeIntervals);
        saccadeTimes = saccadeTimes(...
            saccadeTimes < max(discreteSampleTimes));

        % Find discrete sample time closest to each saccade time
        discreteSaccadeTimes = zeros(size(saccadeTimes));
        saccadeTimeIndices = zeros(size(saccadeTimes));
        for ii = 1:length(saccadeTimes)
            [~, closestIndex] = min(abs(discreteSampleTimes - ...
                saccadeTimes(ii)));
            discreteSaccadeTimes(ii) = discreteSampleTimes(closestIndex);
            saccadeTimeIndices(ii) = closestIndex;
        end 
        % mean(diff(saccadeTimes))
        % mean(diff(discreteSaccadeTimes))

        % Convert to integer sample indices when we'll make a saccade
        %saccadeTimeIndices = round(saccadeTimes / sampTime);

        % Truncate the sample time indices to the duration that
        % we are actually simulating. Need to do this because the list of
        % possible saccade times was created to be very long, at least for
        % reasonable integration times.
        % saccadeTimeIndices = saccadeTimeIndices(saccadeTimeIndices <= ...
        %     nEyeMovements);
        % OLD: saccadeTimeIndices = tPos(1:find(saccadeTimeIndices <= ...
        %     nEyeMovements, 1, 'last'));

        % Make sure we do at least one saccade, and in this case put it in
        % the middle of the time sequence
        %tPos = max(tPos, 1);
        %if (isempty(saccadeTimeIndices))
        %    saccadeTimeIndices = round(nEyeMovements / 2);
        %end

        % Comptue saccades by looping over when they happen
        for ii = 1:length(saccadeTimeIndices)
            % Find out where we are at the time of a saccade
            curPos = pos(saccadeTimeIndices(ii), :);

            % Each saccade has a slightly different speed, compute the
            % speed for this saccade. Speed starts cones/sample. convert to
            % cones/sec for clarity of thinking.
            curSaccadeSpeedConesPerSample = speed + speedSD * randn;
            if (curSaccadeSpeedConesPerSample < 0)
                curSaccadeSpeedConesPerSample = speed;
            end
            curSaccadeSpeedConesPerSec = curSaccadeSpeedConesPerSample ...
                / sampTime;

            % Saccade duration and direction are computed from the current
            % position. The code computes a unit vector headed generally
            % towards the origin, with some angular noise. This code makes
            % our algorithm have the feature that saccade amplitude and
            % duration is determined by the average distance from fixation
            % that drift and tremor produce in a inter-saccade interval.
            % This might be good, but we can't have this feature and
            % control saccade duration and magnitude independently of drift
            % and tremor. Indeed, without drift or tremor, there will be no
            % microsaccades.
            saccadeAngle = atand(curPos(2) / curPos(1)) + dirSD * randn;
            saccadeDirection = [cosd(saccadeAngle) sind(saccadeAngle)];
            saccadeDirection = abs(saccadeDirection) .* ...
                (2 * (curPos < 0) - 1);

            % Compute saccadde duration. This is such that we would end up
            % right back at the origin, if direction were not perturbed.
            % Need to have at least one sample.
            saccadeMagnitude = sqrt(curPos(1) ^ 2 + curPos(2) ^ 2);
            saccadeDurationSecs = saccadeMagnitude / ...
                curSaccadeSpeedConesPerSec;
            saccadeDurationSamples = round(saccadeDurationSecs / sampTime);
            if (saccadeDurationSamples == 0)
                saccadeDurationSamples = 1;
            end

            % Compute the offset we need to add to the positions vector, to
            % put in this saccade. Start with an array of zeros, and
            % compute position offsets over the duration of the saccade.
            offset = zeros(nEyeMovements, 2);
            indx = saccadeTimeIndices(ii):min(saccadeTimeIndices(ii) + ...
                saccadeDurationSamples - 1, nEyeMovements);

            % For each frame of the saccade, add in the offset for that
            % frame. We want to go saccadeMagnitude over
            % saccadeDurationSamples, in the saccade direction.
            offset(indx, 1) = saccadeMagnitude * saccadeDirection(1) / ...
                saccadeDurationSamples;
            offset(indx, 2) = saccadeMagnitude * saccadeDirection(2) / ...
                saccadeDurationSamples;

            % Convert offsets into cumulative saccade position delta, and
            % add into the positions we already had.
            pos = pos + cumsum(offset);
        end
    end

    % Remove eye movements due to drift and tremor
    if (microSaccadesOnlyFlag), pos = pos - driftTremorEMpos; end

    nTrialsPos(nn, :, :) = pos;
    pos = zeros(nEyeMovements, 2);
end

%% Adjustments for the return
% Round to discrete cone steps
% nTrialsPos = round(nTrialsPos);
% 
% % Return in microns as well
% nTrialsPosMicrons = nTrialsPos * obj.patternSampleSize(1) * 1e6;

% The positions in the mosaic is always just 1 sequence of eye movements
%}
