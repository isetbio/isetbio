function [absorptions, photocurrents, LMSfilters, meanCur] = ...
    computeForOISequence(obj, oiSequence, varargin)
% Compute cone absorptions and optionally photocurrents for a @oiSequence
%
% Syntax:
%   [absorptions, photocurrents, LMSfilters] = ...
%        computeForOISequence(obj, oiSequence, varargin)
%
% Description:
%    Compute the cone absorptions and optionally the photocurrents for an
%    @oiSequence. Can also contain LMS filters and the mean current.
%
%    There are several ways to use this function. The simplest is to send
%    in a single oiSequence and a single eye movement sequence.
%
%        coneMosaic.compute(oiSequence)
%
%    It is also possible to run this for a multiple eye movement paths. In
%    that case, the coneMosaic object contains only the data from last eye
%    movement path. The full data set for all the computed eye movements
%    paths are contained in the function returns: absorptions and
%    photocurrent.
%
%        [absorptions, photocurrents] = ...
%             cMosaic.computeForOISequence(oiSequence);
%
%   We control the photon noise by cmosaic.noiseFlag, and the photocurrent
%   noise by cmosaic.os.noiseFlag. These have the options: 'random',
%   'frozen', and 'none'. When 'frozen', you can send in a 'seed'. [Note:
%   May not be fully implemented yet.]
%
% Inputs:
%    obj                - @coneMosaic object
%    oiSequence         - @oiSequence object
%
% Outputs:
%    absorptions        - The cone photon absorptions
%                         (photon counts in integrationTime)
%    photocurrent       - The cone photocurrent
%    LMSfilters         - The LMS Filters
%    meanCur            - The mean current
%
% Optional key/value pairs:
%    'seed'             - Value of random noise seed.(default 1). 
%    'emPaths'          - [N x M x 2] matrix of N eye movement paths, each
%                         with Mx2 eye positions (default empty)
%    'trialBlockSize'   - How many trials each trialBlock should have.
%                         Default: is emtpy, which results in nTrials (no
%                         blocking). This only has an effect with
%                         @coneMosaicHex mosaics and when nTrials>1 and it
%                         is useful with large mosaics x lots of trials, in
%                         which case the absorptions matrix does not fit in
%                         the RAM. If set to -1, the number of trial blocks
%                         is computed automatically based on the number of
%                         cores and system RAM.
%   'interpFilters'     - The L, M, S temporal impulse response functions
%                         to be employed for computing the L, M, and S
%                         photocurrent responses
%   'meanCur'           - The steady state current caused by the mean
%                         absorption rate.
%   'currentFlag'       - Whether to compute photocurrent(default false). 
%   'theExpandedMosaic' - We need an expanded version of the coneMosaic to
%                         deal with eye movements. For multiple calls to
%                         computeForOISequence, we may want to generate it
%                         once and pass it. If it is empty (default), the
%                         expanded version is generated here.
%   'workerID'            If this field is non-empty (default is empty),
%                         the progress of the computation is printed in the
%                         command window along with the workerID (from a
%                         parfor loop).
%   'workDescription'     A string describing the condition
%                         computed. Used only when workerID is non-empty
%                         (default empty).
%   'beVerbose'          - Whether to display infos (default false).
%
% Notes:
%    * TODO: Confirm if implementation referenced in note above is done.
%
% See Also:
%    CONEMOSAIC, COMPUTE, t_simplePhotocurrentComputation

% History:
%    xx/xx/16  NPC  ISETBIO Team 2016
%    02/19/18  jnm  Formatting
%    06/17/18  NPC  Support cone efficiency correction with eccentricity
%
%   Examples:
%{
    % This is an example of how to do this for 1, 000 eye movement paths
    % (trials), each consisting of 100 eye movements.
   nTrials = 1000;
    nEyeMovements = 100;
   emPaths = zeros(instancesBlockSize, nEyeMovements, 2);
   for iTrial = 1:nTrials
        theEMPaths(iTrial , :, :) = cMosaic.emGenSequence(nEyeMovements);
   end
   [absorptions, photocurrents] = cMosaic.computeForOISequence(...
        theOIsequence, 'emPaths', theEMPaths, 'currentFlag', true);

    % The returned absorptions has an extra dimension (the first one) so
    % that we can calculate for multiple eye movement paths. The
    % absorptions from a single eye movement case would be:
   absorptions(thisTrial, :, :, :)

    % The coneMosaic object (obj) has the absorptions from the last trial
    % and dimensions (row, col, time).
%}

%% Parse inputs
p = inputParser;
p.addRequired('oiSequence', @(x)( (isa(x, 'oiSequence')) || (isa(x, 'oiArbitrarySequence')) ));
p.addParameter('emPaths', [], @isnumeric);
p.addParameter('seed', 1, @isnumeric);
p.addParameter('interpFilters', [], @isnumeric);
p.addParameter('meanCur', [], @isnumeric);
p.addParameter('trialBlockSize', [], @isnumeric);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('workerID', [], @isnumeric);
p.addParameter('workDescription', '', @ischar);
p.addParameter('beVerbose', false, @islogical);
p.addParameter('theExpandedMosaic', []);
p.addParameter('stimulusSamplingInterval', [], @isnumeric);
p.parse(oiSequence, varargin{:});

oiSequence  = p.Results.oiSequence;

currentSeed = p.Results.seed;
emPaths     = p.Results.emPaths;
LMSfilters  = p.Results.interpFilters;
meanCur     = p.Results.meanCur;
trialBlockSize    = p.Results.trialBlockSize;
currentFlag       = p.Results.currentFlag;
workerID          = p.Results.workerID;
workDescription   = p.Results.workDescription;
theExpandedMosaic = p.Results.theExpandedMosaic;
beVerbose = p.Results.beVerbose;


if (currentFlag)
    % save noiseFlag
    copyOfConeMosaicNoiseFlag = obj.noiseFlag;
    % compute absorptions without noise
    obj.noiseFlag = 'none';
end

% Set debugTiming to true to examine the timing between oiFrames and
% partial/full absorptions. In this mode, the computation stops and waits
% for the user to hit enter during each step.
debugTiming = false;

oiTimeAxis = oiSequence.timeAxis;
nTimes = numel(oiTimeAxis);
if (oiSequence.length ~= nTimes)
    error('oiTimeAxis and oiSequence must have equal length\n');
end

if (isempty(emPaths))
    if ndims(obj.emPositions) == 3
        emPaths = obj.emPositions;
    elseif ismatrix(obj.emPositions)
        emPaths = reshape(obj.emPositions, [1 size(obj.emPositions)]);
    else 
        error('Bad shape for emPositions data %d\n',size(obj.emPositions));
    end
end

if (isempty(emPaths))
    error(['Either supply an ''emPaths'' key-value pair, or preload ' ...
        'coneMosaic.emPositions']);
end

if (isempty(theExpandedMosaic))
    padRows = max(max(abs(emPaths(:, :, 2))));
    padCols = max(max(abs(emPaths(:, :, 1))));

    % We need a copy of the object because of eye movements. Make it here
    % instead of in coneMosaic.compute(), which is called multiple times.
    obj.absorptions = [];
    obj.current = [];
    if (isa(obj.os, 'osLinear'))
        obj.os.lmsConeFilter = [];
    end
    theExpandedMosaic = obj.copy();
    theExpandedMosaic.pattern = zeros(obj.rows + 2 * padRows, ...
        obj.cols + 2 * padCols);
    
end

% Determine if we need to compute eccentricity-dependent corrections for 
% the absorptions, and if so do it here.  
if (obj.shouldCorrectAbsorptionsWithEccentricity())
if (isempty(obj.coneEfficiencyCorrectionFactors))
    correctionFactors = coneMosaicHex.computeConeEfficiencyCorrectionFactors(obj, ...
        mfilename(), 'beVerbose', beVerbose);
    obj.setConeQuantalEfficiencyCorrectionFactors(correctionFactors);
end
end

%% Get ready for output variables
photocurrents = [];
% varargout{1} = [];
% varargout{2} = [];

%% Save default integration time
defaultIntegrationTime = obj.integrationTime;

%% Compute eye movement time axis
nTrials = size(emPaths, 1);
nEyeMovements = size(emPaths, 2);

%% round everything to 1 microsecond to avoid rounding-off timing issues
rounded.factor = 1 / (1000 * 1000.0);
rounded.oiTimeAxis = round(oiTimeAxis / rounded.factor);
rounded.defaultIntegrationTime = round(defaultIntegrationTime / ...
    rounded.factor);
rounded.eyeMovementTimeAxis = rounded.oiTimeAxis(1) + ...
    (0:1:(nEyeMovements - 1)) * rounded.defaultIntegrationTime;

%% Compute OIrefresh
if (numel(rounded.oiTimeAxis) == 1)
    if (~isempty(p.Results.stimulusSamplingInterval))
        rounded.oiRefreshInterval = ...
           round(p.Results.stimulusSamplingInterval/rounded.factor);
    else
        % No information about what the stimulus refresh interval is,
        % so arbitrarily set it to the integrationTime
        rounded.oiRefreshInterval = rounded.defaultIntegrationTime;
    end
else
    rounded.oiRefreshInterval = rounded.oiTimeAxis(2) - ...
        rounded.oiTimeAxis(1);
end

% Only allocate memory for the non-null cones in a 3D matrix [instances x
% numel(nonNullConesIndices) x time]
nonNullConesIndices = find(obj.pattern > 1);
absorptions = zeros(nTrials, numel(nonNullConesIndices), ...
    numel(rounded.eyeMovementTimeAxis), 'single');

% Organize trials in blocks if we have a hex mosaic
if (isempty(trialBlockSize)), trialBlockSize = nTrials; end
if (nTrials > 1) && (trialBlockSize == -1) && (isa(obj, 'coneMosaicHex'))
    blockedTrialIndices = ...
        computeBlockedTrialIndices(trialBlockSize, nTrials);
else
    blockedTrialIndices{1} = 1:nTrials;
end

if (rounded.oiRefreshInterval >= rounded.defaultIntegrationTime)
    % There are two main time sampling scenarios. This one is when the oi
    % update rate is SLOWER than the cone integration time which is also
    % equal to the eye movement update rate.
    %
    
    %% SCENARIO
    % oiRI = oiRefreshInterval
    % p1, p2 = partial absorption
    % full = full absorption
    %
    %          | ----- oiRI ----- | ----- oiRI ----- | ----- oiRI ----- |
    %    oi    |                  |                  |                  |
    % sequence |                  |<---tFrameStart   |                  |
    %__________.__________________|     tFrameEnd--->|__________________|
    %   eye    |     |     |        |        |        |     |     |     |
    % movement |     |     |intTime>|intTime>|intTime>|     |     |     |
    % sequence |     |     |        |        |        |     |     |     |
    %----------------------|****xxxx|********|********|-----------------|
    %                       p1    p2   full     full

    if (debugTiming)
        hFig = figure(9000);
        clf;
        set(hFig, 'Name', 'Visual Timing Debug Window', ...
            'Position', [10 10 1600 200]);
        subplot('Position', [0.01 0.01 0.98 0.98]);
        hold on
    end

    % Loop over the optical images
    for oiIndex = 1:oiSequence.length
        if (~isempty(workerID))
            % Update progress in command window
            displayProgress(workerID, ...
                sprintf('%s-absorptions', workDescription), ...
                0.5 * oiIndex / oiSequence.length);
        end

        % Precompute the two OIs to save computation time
        previousOI = [];
        if (oiIndex > 1)
            previousOI = oiSequence.frameAtIndex(oiIndex - 1);
        end
        currentOI = oiSequence.frameAtIndex(oiIndex);

        % Current oi time limits
        tFrameStart = rounded.oiTimeAxis(oiIndex);
        tFrameEnd = tFrameStart + rounded.oiRefreshInterval;

        if (debugTiming)
            x = [tFrameStart tFrameStart tFrameEnd tFrameEnd tFrameStart];
            y = [0 -10 -10 0 0];
            plot(x, y, 'r-');
        end

        % Find eye movement indices withing the oi limits
         indices = find( ...
            (rounded.eyeMovementTimeAxis > tFrameStart - ...
            rounded.defaultIntegrationTime) & ...
            (rounded.eyeMovementTimeAxis <= tFrameEnd  - ...
            rounded.defaultIntegrationTime));

        if (isempty(indices))
            if (debugTiming)
                fprintf('No eye movements within oiIndex #%d\n', oiIndex);
            end
            continue;
        end

        for iTrialBlock = 1:numel(blockedTrialIndices)
            % Obtain trials to be computed at this trial block
            trialIndicesForBlock = blockedTrialIndices{iTrialBlock};

            % the first eye movement requires special treatment as it may
            % have started before the current frame, so we need to compute
            % partial absorptions over the previous frame and over the
            % current frame
            idx = indices(1);
            integrationTimeForFirstPartialAbsorption = tFrameStart - ...
                rounded.eyeMovementTimeAxis(idx);
            integrationTimeForSecondPartialAbsorption = ...
                rounded.eyeMovementTimeAxis(idx) + ...
                rounded.defaultIntegrationTime - tFrameStart;
            
            if (integrationTimeForFirstPartialAbsorption < ...
                    eps(tFrameStart))
                integrationTimeForFirstPartialAbsorption = 0;
            elseif (integrationTimeForSecondPartialAbsorption < ...
                    eps(rounded.eyeMovementTimeAxis(idx) + ...
                    rounded.defaultIntegrationTime))
                integrationTimeForSecondPartialAbsorption = 0;
            end

            % Partial absorptions (p1 in graph above) with previous oi
            % (across all instances)
            if (oiIndex > 1) && ...
                    (integrationTimeForFirstPartialAbsorption > 0)
                % Update the @coneMosaic with the partial integration time
                obj.integrationTime = ...
                    integrationTimeForFirstPartialAbsorption * ...
                    rounded.factor;
                % Compute partial absorptions
                emSubPath = reshape(squeeze(emPaths(...
                    trialIndicesForBlock, idx, :)), ...
                    [numel(trialIndicesForBlock) 2]);
                currentSeed = currentSeed  + 1;
                absorptionsAllTrials = single(obj.compute(...
                    previousOI, ...
                    'theExpandedMosaic', theExpandedMosaic, ...
                    'seed', currentSeed , ...
                    'emPath', emSubPath, ...
                    'currentFlag', false, ...
                    'beVerbose', beVerbose));
                if (debugTiming)
                    x = tFrameStart - [0 0 obj.integrationTime / ...
                        rounded.factor obj.integrationTime / ...
                        rounded.factor 0];
                    y = [-1 -2 -2 -1 -1] - (oiIndex - 1) * 0.1;
                    plot(x, y, 'b-', 'LineWidth', 1.5);
                end
            else
                absorptionsAllTrials = zeros(size(obj.pattern, 1), ...
                    size(obj.pattern, 2), numel(trialIndicesForBlock), ...
                    'single');
            end

            if (integrationTimeForSecondPartialAbsorption > 0)
                % Partial absorptions (p2 in graph above) with current oi
                % across all instances. Update the @coneMosaic with the
                % partial integration time
                obj.integrationTime = ...
                    integrationTimeForSecondPartialAbsorption * ...
                    rounded.factor;
                % Compute partial absorptions
                emSubPath = reshape(squeeze(...
                    emPaths(trialIndicesForBlock, idx, :)), ...
                    [numel(trialIndicesForBlock) 2]);
                currentSeed = currentSeed  + 1;
                absorptionsAllTrials =  absorptionsAllTrials + ...
                    single(obj.compute(currentOI, ...
                    'theExpandedMosaic', theExpandedMosaic, ...
                    'seed', currentSeed, ...
                    'emPath', emSubPath, ...
                    'currentFlag', false, ...
                    'beVerbose', beVerbose));
                if (debugTiming)
                    x = tFrameStart + [0 0 obj.integrationTime / ...
                        rounded.factor obj.integrationTime / ...
                        rounded.factor 0];
                    y = [-1 -2 -2 -1 -1] - (oiIndex - 1) * 0.1;
                    plot(x, y, 'm-', 'LineWidth', 1.5);
                end
                % vcNewGraphWin;
                % imagesc(absorptionsDuringCurrentFrame);
            end

            % Reformat and insert to time series
            firstEMinsertionIndex = round((...
                rounded.eyeMovementTimeAxis(idx) - ...
            rounded.eyeMovementTimeAxis(1)) / ...
            rounded.defaultIntegrationTime) + 1;
            reformatAbsorptionsAllTrialsMatrix(numel(...
                trialIndicesForBlock), 1, size(obj.pattern, 1), ...
                size(obj.pattern, 2));
            absorptions(trialIndicesForBlock, :, ...
                firstEMinsertionIndex) = absorptionsAllTrials;
            lastInsertionIndex = firstEMinsertionIndex;

            % Full absorptions with current oi and default integration time
            if (numel(indices)>1)
                % Update the @coneMosaic with the default integration time
                obj.integrationTime = rounded.defaultIntegrationTime * ...
                    rounded.factor;
                % Compute absorptions for all remaining the eye movements
                idx = indices(2:end);
                emSubPath = reshape(...
                    emPaths(trialIndicesForBlock, idx, :), ...
                    [numel(trialIndicesForBlock) * numel(idx) 2]);
                currentSeed = currentSeed + 1;
                absorptionsAllTrials = single(obj.compute(...
                    currentOI, ...
                    'theExpandedMosaic', theExpandedMosaic, ...
                    'seed', currentSeed, ...
                    'emPath', emSubPath, ...
                    'currentFlag', false, ...
                    'beVerbose', beVerbose));
                % Reformat and insert to time series
                remainingEMinsertionIndices = round((...
                    rounded.eyeMovementTimeAxis(idx) - ...
                    rounded.eyeMovementTimeAxis(1)) / ...
                    rounded.defaultIntegrationTime) + 1;
                reformatAbsorptionsAllTrialsMatrix(...
                    numel(trialIndicesForBlock), ...
                    numel(remainingEMinsertionIndices), ...
                    size(obj.pattern, 1), size(obj.pattern, 2));
                absorptions(trialIndicesForBlock, :, ...
                    remainingEMinsertionIndices) = absorptionsAllTrials;
                lastInsertionIndex = remainingEMinsertionIndices(end);
                
                if (debugTiming)
                    for kk = 2:numel(indices)
                        x = tFrameStart + ...
                            integrationTimeForSecondPartialAbsorption + ...
                            [0, 0, (kk - 1) * obj.integrationTime / ...
                            rounded.factor, (kk - 1) * ...
                            obj.integrationTime / rounded.factor, 0];
                        y = [-1 -2 -2 -1 -1] - (oiIndex - 1) * 0.1;
                        plot(x, y, 'k-', 'LineWidth', 1.5);
                        drawnow;
                    end
                end
            end
        end % iTrialBlock

        % If we are in the last OIframe, and there are still some
        % EmInsertionIndices that are unfilled, fill them with absorptions
        % from the last computed insertion index
        if (oiIndex == oiSequence.length) && ...
                (lastInsertionIndex < numel(rounded.eyeMovementTimeAxis))
           if (lastInsertionIndex == ...
                   numel(rounded.eyeMovementTimeAxis) - 1)
            % fill last entry
            absorptions(trialIndicesForBlock, :, ...
                lastInsertionIndex + 1:end) = ...
            absorptions(trialIndicesForBlock, :, lastInsertionIndex);
           else
               % Return an error if more than one emIndices are nonfilled,
               % as this should not be hapening.
               error(['oiIndex: %d. empty absorption matrix indices to' ...
                   ' be filled: %d\n'], oiIndex, ...
                   size(absorptions, 3) - lastInsertionIndex);
           end
        end
        
        if (debugTiming)
            if(numel(indices) == 1)
                fprintf(['[%3d/%3d]: p1=%05.1fms p2=%05.1fms ' ...
                    'f=%05.1fms, i1=%03d, iR=[]      i1Time=%06.1fms, ' ...
                    'iRtime=[]                     (%d)\n'], oiIndex, ...
                    oiSequence.length, ...
                    integrationTimeForFirstPartialAbsorption * ...
                    rounded.factor * 1000, ...
                    integrationTimeForSecondPartialAbsorption * ...
                    rounded.factor * 1000, 0, firstEMinsertionIndex, ...
                    rounded.eyeMovementTimeAxis(indices(1)) * ...
                    rounded.factor * 1000, numel(indices));
            else
                fprintf(['[%3d/%3d]: p1=%05.1fms p2=%05.1fms ' ...
                    'f=%05.1fms, i1=%03d, iR=%03d-%03d i1Time=%06.1fms' ...
                    ', iRtime=[%06.1fms .. %06.1fms] (%d)\n'], oiIndex, ...
                    oiSequence.length, ...
                    integrationTimeForFirstPartialAbsorption * ...
                    rounded.factor * 1000, ...
                    integrationTimeForSecondPartialAbsorption * ...
                    rounded.factor * 1000, obj.integrationTime * 1000, ...
                    firstEMinsertionIndex, ...
                    remainingEMinsertionIndices(1), ...
                    remainingEMinsertionIndices(end), ...
                    rounded.eyeMovementTimeAxis(indices(1)) * ...
                    rounded.factor * 1000, ...
                    rounded.eyeMovementTimeAxis(indices(2)) * ...
                    rounded.factor * 1000, ...
                    rounded.eyeMovementTimeAxis(indices(end)) * ...
                    rounded.factor * 10000, numel(indices));
            end
            if (oiSequence.length > 1)
                set(gca, 'XLim', [rounded.oiTimeAxis(1), ...
                    rounded.oiTimeAxis(end) + ...
                    integrationTimeForSecondPartialAbsorption]);
                pause
            end
        end
    end  % oiIndex

    % Clear OIs - not needed anymore
    varsToClear = {'currentOI', 'previousOI'};
    clear(varsToClear{:});

    % rounded.oiRefreshInterval > rounded.defaultIntegrationTime
else
    % There are two main time sampling scenarios. This one is when the oi
    % update rate is FASTER than the cone integration time which is also
    % equal to the eye movement update rate
    
    %% SCENARIO
    % oiR = oiRefresh
    % p1, p2 = partial absorption
    % full = full absorption
    %
    %   eye    |               |                        |               |
    % movement |               |                        |               |
    % sequence |               |<----emStart            |               |
    %__________._______________|            emEnd------>|_______________|
    %    oi    |     |     |---oiR--->|---oiR--->|---oiR--->|     |     |
    % sequence |     |     |    ******|**********|******    |     |     |
    %-------------------------------p1-----full------p2------------------

    % Loop over the eye movements
    for emIndex = 1:nEyeMovements
        if (~isempty(workerID))
            displayProgress(workerID, ...
                sprintf('%s-absorptions', workDescription), ...
                0.5 * emIndex / nEyeMovements);
        end

        % Current eye movement time limits
        emStart = rounded.eyeMovementTimeAxis(emIndex);
        emEnd = emStart + rounded.defaultIntegrationTime;

        % Find oi indices withing the eye movement frame time limits
        indices = find(( ...
            rounded.oiTimeAxis > emStart - rounded.oiRefreshInterval) & ...
            (rounded.oiTimeAxis <= emEnd));

        if (isempty(indices))
            if (debugTiming)
                fprintf('No OIs within emIndex #%d\n', emIndex);
            end
            continue;
        end

        for iTrialBlock = 1:numel(blockedTrialIndices)
            % Obtain trials to be computed at this trial block
            trialIndicesForBlock = blockedTrialIndices{iTrialBlock};
            
            % sum the partial integration times
            actualIntegrationTime = 0;
            
            % Partial absorptions during the ovelap with the OI that
            % started before the emStart
            idx = indices(1);
            integrationTimeForFirstPartialAbsorption = ...
                rounded.oiTimeAxis(idx) + rounded.oiRefreshInterval - ...
                emStart;
            if (integrationTimeForFirstPartialAbsorption > ...
                    eps(rounded.oiTimeAxis(idx) + ...
                    rounded.oiRefreshInterval))
                % Update the @coneMosaic with the partial integration time
                obj.integrationTime = ...
                    integrationTimeForFirstPartialAbsorption * ...
                    rounded.factor;
                % Update the sum of partial integration times
                actualIntegrationTime = actualIntegrationTime + ...
                    obj.integrationTime;
                % Compute absorptions
                emSubPath = reshape(emPaths(trialIndicesForBlock, ...
                    emIndex, :), [numel(trialIndicesForBlock) 2]);
                currentSeed = currentSeed  + 1;
                absorptionsAllTrials = single(obj.compute(...
                    oiSequence.frameAtIndex(idx), ...
                    'theExpandedMosaic', theExpandedMosaic, ...
                    'seed', currentSeed, ...
                    'emPath', emSubPath, ...
                    'currentFlag', false));
            else
                absorptionsAllTrials = zeros(size(obj.pattern, 1), ...
                    size(obj.pattern, 2), numel(trialIndicesForBlock), ...
                    'single');
            end

            % Next, compute full absorptions for the remaining OIs
            if (numel(indices) > 1)
                for k = 2:numel(indices)
                    % Update @coneMosaic with the partial integration time
                    if (k < numel(indices))
                        % full absorption during the duration of the OI
                        obj.integrationTime = rounded.oiRefreshInterval ...
                            * rounded.factor;
                    else
                        % partial absorption during the last OI
                        idx = indices(end);
                        integrationTimeForLastPartialAbsorption = emEnd ...
                            - rounded.oiTimeAxis(idx);
                        obj.integrationTime = ...
                            integrationTimeForLastPartialAbsorption * ...
                            rounded.factor;
                    end

                    % Update the sum of partial integration times
                    actualIntegrationTime = actualIntegrationTime + ...
                        obj.integrationTime;
                    % Compute absorptions
                    emSubPath = reshape(emPaths(trialIndicesForBlock, ...
                        emIndex, :), [numel(trialIndicesForBlock) 2]);
                    currentSeed = currentSeed  + 1;
                    absorptionsAllTrials = absorptionsAllTrials + ...
                        single(obj.compute(oiSequence.frameAtIndex(...
                        indices(k)), 'theExpandedMosaic', ...
                        theExpandedMosaic, 'seed', currentSeed, ...
                        'emPath', emSubPath, ...
                        'currentFlag', false));
                end  % for k
            end

            % Check the sum of partial absorptions equal the default
            % absorption time.
            if (abs(actualIntegrationTime - ...
                    rounded.defaultIntegrationTime * rounded.factor) > ...
                    rounded.defaultIntegrationTime * 0.0001) && ...
                    (emIndex < numel(rounded.eyeMovementTimeAxis))
                error(['Actual integration time (%3.5f) not equal to ' ...
                    'desired value (%3.5f) [emIndex: %d / %d]\n'], ...
                    actualIntegrationTime, ...
                    rounded.defaultIntegrationTime, emIndex, ...
                    numel(rounded.eyeMovementTimeAxis));
            end

            % Reformat and insert to time series
            insertionIndices = round((...
                rounded.eyeMovementTimeAxis(emIndex) - ...
                rounded.eyeMovementTimeAxis(1)) / ...
                rounded.defaultIntegrationTime) + 1;
            reformatAbsorptionsAllTrialsMatrix(...
                numel(trialIndicesForBlock), numel(insertionIndices), ...
                size(obj.pattern, 1), size(obj.pattern, 2));
            absorptions(trialIndicesForBlock, :, insertionIndices) = ...
                absorptionsAllTrials;
        end % iTrialBlock

    end % emIndex
end % rounded.oiRefreshInterval > rounded.defaultIntegrationTime

% Save some RAM
clear 'oiSequence';

% Restore default integrationTime
obj.integrationTime = defaultIntegrationTime;

% Reload the full eye movement sequence for the last trial only
if (ndims(emPaths) == 3)
    obj.emPositions = reshape(squeeze(emPaths(nTrials, :, :)), ...
        [nEyeMovements 2]);
else
    obj.emPositions = emPaths;
end

% If we don't compute the current, we return only the absorptions from the
% last trial. We never compute the current if there are no eye movements.
if (~currentFlag) || (numel(rounded.eyeMovementTimeAxis) == 1)
    returnAbsorptionsFromLastTrial();
    return;
end

%% Photocurrent computation
% The currentFlag must be on, and there must be a few eye movements. So we
% compute. N.B. Not adequately checked for osBioPhys model. Runs ok for
% osLinear model.

% free some RAM before procedding with allocating more RAM
obj.absorptions = [];

if (isa(obj, 'coneMosaicHex'))
    photocurrents = zeros(nTrials, numel(nonNullConesIndices), ...
        numel(rounded.eyeMovementTimeAxis), 'single');
    for ii = 1:nTrials
        if (~isempty(workerID)) && (mod(ii - 1, round(nTrials / 10)) == 0)
            displayProgress(workerID, sprintf('%s-current', ...
                workDescription), 0.5 + 0.5 * ii / nTrials);
        end

        currentSeed = currentSeed  + 1;
        if ii == 1 && (isempty(meanCur) || isempty(LMSfilters))
            % On the first trial, compute the interpolated linear
            % filters and the mean current, unless they were passed in
            [LMSfilters, meanCur] = obj.computeCurrent(...
                'seed', currentSeed, ...
                'absorptionsInXWFormat', squeeze(absorptions(ii, :, :)));
        else
            LMSfilters = obj.computeCurrent(...
                'seed', currentSeed, ...
                'interpFilters', LMSfilters, ...
                'meanCur', meanCur, ...
                'absorptionsInXWFormat', squeeze(absorptions(ii, :, :)));
        end

        photocurrents(ii, :, :) = single(obj.current);
    end
else
    photocurrents = zeros(nTrials, obj.rows, obj.cols, ...
        numel(rounded.eyeMovementTimeAxis), 'single');
    for ii = 1:nTrials
        if (~isempty(workerID)) && (mod(ii, round(nTrials / 10)) == 0)
            displayProgress(workerID, sprintf('%s-current', ...
                workDescription), 0.5 + 0.5 * ii / nTrials);
        end
        % Put this trial of absorptions into the cone mosaic
        obj.absorptions = reshape(squeeze(absorptions(ii, :, :, :)), ...
            [obj.rows obj.cols, numel(rounded.eyeMovementTimeAxis)]);
        currentSeed = currentSeed  + 1;
        if ii == 1 && (isempty(meanCur) || isempty(LMSfilters))
            % On the first trial, compute the interpolated linear
            % filters and the mean current
            [LMSfilters, meanCur] = obj.computeCurrent('seed', ...
                currentSeed);
        else
            % On the remaining trials, use the same interpolated
            % filters and mean current
            LMSfilters = obj.computeCurrent('seed', currentSeed, ...
                'interpFilters', LMSfilters, 'meanCur', meanCur);
        end
        photocurrents(ii, :, :, :) = single(reshape(obj.current, ...
            [1 obj.rows obj.cols numel(rounded.eyeMovementTimeAxis)]));
    end
end

%% Since we computed noise-free absorptions (so that the photocurrent
%% is computed on the noise-free absorptions), we need to add photon noise
%% to the absorptions at this point

obj.noiseFlag = copyOfConeMosaicNoiseFlag;
if ~(strcmp(obj.noiseFlag, 'none'))
    %% add photon noise
    absorptions = obj.photonNoise(absorptions, ...
         'noiseFlag', obj.noiseFlag, 'seed', currentSeed);
end



% Reload the absorptions from the last instance (again, since we destroyed
% obj.absorptions in the current computation)
returnAbsorptionsFromLastTrial();

%% Nested function for returning absorptions from the last trial
    function returnAbsorptionsFromLastTrial()
        % Return the absorptions from the last trial
        %
        % Syntax:
        %   returnAbsorptionsFromLastTrial()
        %
        % Definition:
        %    Return the absorptions from the last trial.
        %
        % Inputs:
        %    None.
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    None.
        %
        if (isa(obj, 'coneMosaicHex'))
            tmp = squeeze(absorptions(nTrials, :, :));
            if (numel(rounded.eyeMovementTimeAxis) == 1), tmp = tmp'; end
            
            % Return the absorptions from the last trial after reshaping to
            % full 3D matrix [cone_rows, cone_cols, time]
            obj.absorptions = obj.reshapeHex2DmapToHex3Dmap(tmp);
        else
            % Reshape to full 4D matrix
            % [instances, cone_rows, cone_cols, time]
            absorptions = reshape(absorptions, [nTrials, ...
                size(obj.pattern, 1), size(obj.pattern, 2), ...
                numel(rounded.eyeMovementTimeAxis)]);
            
            % Return the absorptions from the last trial.
            obj.absorptions = reshape(...
                squeeze(absorptions(nTrials, :, :)), ...
                [size(obj.pattern, 1), size(obj.pattern, 2), ...
                numel(rounded.eyeMovementTimeAxis)]);
        end
    end

%% Nested function to reformat absorptions
    function reformatAbsorptionsAllTrialsMatrix(trialsNum, ...
            timePointsNum, coneRows, coneCols)
        % Reformat all of the trial matrix absorptions
        %
        % Syntax:
        %   reformatAbsorptionsAllTrialsMatrix(trialsNum, ...
        %       timePointsNum, coneRows, coneCols)
        %
        % Description:
        %    This is a nested function to reformat the absorptions
        %
        % Inputs:
        %    trialsNum     - The number of trials
        %    timePointsNum - The number of time points
        %    coneRows      - The cone rows
        %    coneCols      - The cone columns
        %
        % Outputs:
        %    None
        %
        % Optional key/value pairs:
        %    None.
        %

        % Reshape to cones x instances x timePoints. Note the 3rd dimension
        % of absorptionsAllTrials is traditionally time, but when we are
        % doing multiple instances, it is instances * time
        absorptionsAllTrials = reshape(absorptionsAllTrials, ...
            [coneRows * coneCols, trialsNum, timePointsNum]);

        % Only get the absorptions for the non-null cones - this has an
        % effect only for coneMosaicHex mosaics
        absorptionsAllTrials = absorptionsAllTrials(...
            nonNullConesIndices, :, :);

        % Reshape to [instances x cones x timePoints]
        absorptionsAllTrials = permute(absorptionsAllTrials, [2 1 3]);
    end
end

%%
function displayProgress(workerID, workDescription, progress)
% Display the progress
%
% Syntax:
%   displayProgress(workerID, workDescription, progress)
%
% Description:
%    Display the work progress
%
% Inputs:
%    workerID - Integer. Worker ID
%    workDescription - String. Describe the work
%    progress - Double. Percentage complete.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None
%
displayFormat = 'numeric';
if (strcmp(displayFormat, 'numeric'))
    if (isnan(progress))
        fprintf('worker-%02d: %s \n', workerID, workDescription);
    else
        fprintf('worker-%02d: %s (%2.0f%%)\n', workerID, ...
            workDescription, progress * 100);
    end
else
    maxStarsNum = 32;
    if (isnan(progress))
        fprintf('worker-%02d: %s |', workerID, workDescription);
        for k = 1:maxStarsNum, fprintf('*'); end
        fprintf('|');
    else
        fprintf('worker-%02d: %s |', workerID, workDescription);
        if (progress > 1), progress = 1; end
        starsNum = round(maxStarsNum * progress);
        for k = 1:starsNum, fprintf('*'); end
    end
    fprintf('\n');
end

end