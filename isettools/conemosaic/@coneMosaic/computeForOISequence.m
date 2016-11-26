function [absorptions, photocurrents] = computeForOISequence(obj, oiSequence, varargin)
% Compute cone absorptions and (optionally) photocurrents for a @oiSequence
%
%    [absorptions, photocurrents] = cMosaic.compute(oiSequence, varargin);
%
% Inputs:
%   obj         - @coneMosaic object
%   oiSequence  - @oiSequence object
%
% Optional inputs:
%   emPaths      - [N x M x 2] matrix of N eye movement paths, each with
%                  Mx2 eye positions
%   currentFlag  - logical, whether to compute photocurrent
%   newNoise     - logical, whether to use new random seed
%
% Outputs:
%   absorptions          - cone photon absorptions (photon counts in integrationTime)
%   photocurrent         - cone photocurrent
%
% There are several ways to use this function.  The simplest is to send in
% a single oiSequence and a single eye movement sequence.
%
%   coneMosaic.compute(oiSequence)
%
% It is also possible to run this for a multiple eye movement paths. In
% that case, the coneMosaic object contains only the last trial.  The full
% set of data for all the trial are contained in the returned outputs,
% absorptions and photocurrent.
%
%   [absorptions, photocurrents] = cMosaic.computeForOISequence(oiSequence);
%
% Examples:
%  This is an example of how to do this for 1,000 eye movement paths
%  (trials), each consisting of 100 eye movements.
%
%   nTrials = 1000; eyeMovementsNum = 100; 
%   emPaths = zeros(instancesBlockSize, eyeMovementsNum, 2);
%   for iTrial = 1:nTrials
%    theEMPaths(iTrial , :,:) = cMosaic.emGenSequence(eyeMovementsNum);
%   end
%   [absorptions, photocurrents] = cMosaic.computeForOISequence(...
%       theOIsequence, ...
%       'emPaths', theEMPaths, ...
%       'currentFlag', true, ...
%       'newNoise', true ...
%   );
%                    
% The returned absorptions has an extra dimension (the first one) so that
% we can calculate for multiple eye movement paths.  The absorptions from a
% single eye movement case would be
%
%    absorptions(thisTrial,:,:,:)
%
% The coneMosaic object (obj) has the absorptions from the last trial and
% dimensions (row,col,time).
%
% See also: coneMosaic.compute, v_cmosaic, 
%
% NPC ISETBIO Team 2016

%% Parse inputs
p = inputParser;
p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
p.addParameter('emPaths', [], @isnumeric);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('newNoise', true, @islogical);
p.addParameter('workerID', [], @isnumeric);
p.addParameter('workDescription', '', @ischar);
p.parse(oiSequence, varargin{:});

oiSequence  = p.Results.oiSequence;
emPaths     = p.Results.emPaths;
currentFlag = p.Results.currentFlag;
newNoise    = p.Results.newNoise;
workerID    = p.Results.workerID;
workDescription = p.Results.workDescription;
oiTimeAxis  = oiSequence.timeAxis;

nTimes = numel(oiTimeAxis);
if (oiSequence.length ~= nTimes)
    error('oiTimeAxis and oiSequence must have equal length\n');
end

if (isempty(emPaths))
    emPaths = reshape(obj.emPositions, [1 size(obj.emPositions)]);
end

if (isempty(emPaths))
    error('Either supply an ''emPaths'' key-value pair, or preload coneMosaic.emPositions');
end

%% Get ready for output variables

photocurrents = [];
% varargout{1} = [];
% varargout{2} = [];

%% Save default integration time
defaultIntegrationTime = obj.integrationTime;

%% Compute eye movement time axis
nTrials         = size(emPaths,1);
eyeMovementsNum = size(emPaths,2);
eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(eyeMovementsNum-1)) * obj.integrationTime;

%% Compute OIrefresh
if (numel(oiTimeAxis) == 1)
    oiRefreshInterval = defaultIntegrationTime;
else
    oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);
end


% Only allocate memory for the non-null cones in a 3D matrix [instances x numel(nonNullConesIndices) x time]
nonNullConesIndices = find(obj.pattern>1);
absorptions = zeros(nTrials, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');

if (oiRefreshInterval >= defaultIntegrationTime)
    % There are two main time sampling scenarios.  This one is when the oi
    % update rate is SLOWER than the cone integration time which is also
    % equal to the eye movement update rate.  
    % 
    
    % SCENARIO
    %              |----- oiRefreshInterval ----|----- oiRefreshInterval ----|----- oiRefreshInterval ----|
    %      oi      |                            |                            |                            |
    %   sequence   |                            |<-----tFrameStart           |                            |
    % _____________.____________________________|           tFrameEnd------->|____________________________|
    %     eye      |          |          |          |          |          |          |          |          |
    %   movement   |          |          |-intTime->|-intTime->|-intTime->|          |          |          |
    %   sequence   |          |          |          |          |          |          |          |          |
    % ------------------------------------******xxxx|**********|**********|--------------------------------|
    %                                       p1   p2     full        full
    %                   partial absorption_/      \_partial absorption  \_ full absorption
    %
    
    % Loop over the optical images
    for oiIndex = 1:oiSequence.length
        
        if (~isempty(workerID))
            % Update progress in command window
            displayProgress(workerID, workDescription, 0.5*oiIndex/oiSequence.length);
        end
        
        % Current oi time limits
        tFrameStart = oiTimeAxis(oiIndex);
        tFrameEnd   = tFrameStart + oiRefreshInterval;
        
        % Find eye movement indices withing the oi limits
        indices = find( (eyeMovementTimeAxis >=  tFrameStart-defaultIntegrationTime) & ...
            (eyeMovementTimeAxis <= tFrameEnd-defaultIntegrationTime+eps) );
        
        if (isempty(indices))
            % time samples in
            % the mosaic than we have oi samples.  That should be OK.
            % disp('Fewer Eye movement time samples than oi samples')
            % break;
            error('Empty indices. This should never happen.');
        end
        % the first eye movement requires special treatment as it may have started before the current frame,
        % so we need to compute partial absorptions over the previous frame and over the current frame
        idx = indices(1);
        integrationTimeForFirstPartialAbsorption = tFrameStart-eyeMovementTimeAxis(idx);
        integrationTimeForSecondPartialAbsorption = eyeMovementTimeAxis(idx)+defaultIntegrationTime-tFrameStart;
        
        % Partial absorptions (p1 in graph above) with previous oi
        % (across all instances)
        if (oiIndex > 1)
            % Update the @coneMosaic with the partial integration time
            obj.integrationTime = integrationTimeForFirstPartialAbsorption;
            % Compute partial absorptions
            emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
            obj.absorptions = [];
            absorptionsDuringPreviousFrame = obj.compute(...
                oiSequence.frameAtIndex(oiIndex-1), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
        else
            absorptionsDuringPreviousFrame = zeros(size(obj.pattern,1), size(obj.pattern,2), nTrials);
        end
        
        % Partial absorptions (p2 in graph above) with current oi across
        % all instances.
        % Update the @coneMosaic with the partial integration time
        obj.integrationTime = integrationTimeForSecondPartialAbsorption;
        
        % Compute partial absorptions
        emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
        obj.absorptions = [];
        absorptionsDuringCurrentFrame =  obj.compute(...
            oiSequence.frameAtIndex(oiIndex), ...
            'emPath', emSubPath, ...
            'newNoise', newNoise, ...
            'currentFlag', false ...
            );
        % vcNewGraphWin; imagesc(absorptionsDuringCurrentFrame);
        
        % summed absorptions
        absorptionsAllTrials = ...
            absorptionsDuringPreviousFrame + absorptionsDuringCurrentFrame;
        
        % Reformat and insert to time series
        insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
        absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
        
        % Full absorptions with current oi and default integration time
        if (numel(indices)>1)
            % Update the @coneMosaic with the default integration time
            obj.integrationTime = defaultIntegrationTime;
            % Compute absorptions for all remaining the OIs
            idx = indices(2:end);
            emSubPath = reshape(emPaths(1:nTrials, idx,:), [nTrials*numel(idx) 2]);
            obj.absorptions = [];
            absorptionsAllTrials = obj.compute(...
                oiSequence.frameAtIndex(oiIndex), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
            % Reformat and insert to time series
            insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
            absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
        end
    end  % oiIndex
    
    % oiRefreshInterval > defaultIntegrationTime
else
    % There are two main time sampling scenarios.  This one is when the oi
    % update rate is FASTER than the cone integration time which is also equal
    % to the eye movement update rate
    
    % SCENARIO
    %     eye      |                             |                             |                             |
    %   movement   |                             |                             |                             |
    %   sequence   |                             |<-----emStart                |                             |
    % _____________._____________________________|                emEnd------->|_____________________________|
    %              |            |            |            |            |            |          |          |
    %      oi      |            |            |-oiRefresh->|-oiRefresh->|-oiRefresh->|          |          |
    %   sequence   |            |            |    ********|************|*******     |          |          |
    % ---------------------------------------------partial-----full-----partial-------------------------------
    %                                            absorption  absorption absorption
    
    % Loop over the eye movements
    for emIndex = 1:eyeMovementsNum
        
        if (~isempty(workerID))
            displayProgress(workerID, workDescription, 0.5*emIndex/eyeMovementsNum);
        end
        
        % Current eye movement time limits
        emStart = eyeMovementTimeAxis(emIndex);
        emEnd   = emStart + defaultIntegrationTime;
        
        % sum the partial integration times
        actualIntegrationTime = 0;
        
        % Find oi indices withing the eye movement frame time limits
        indices = find( (oiTimeAxis >= emStart-oiRefreshInterval) & ...
            (oiTimeAxis <= emEnd + eps) );
        
        if (isempty(indices))
            %disp('Fewer eye movement time samples than oi samples')
            %break;
            error('Empty indices. This should never happen.');
        end
        
        % Partial absorptions during the ovelap with the OI that started before the emStart
        idx = indices(1);
        integrationTimeForFirstPartialAbsorption = oiTimeAxis(idx)+oiRefreshInterval-emStart;
        % Update the @coneMosaic with the partial integration time
        obj.integrationTime = integrationTimeForFirstPartialAbsorption;
        % Update the sum of partial integration times
        actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
        % Compute absorptions
        emSubPath = reshape(emPaths(1:nTrials, emIndex,:), [nTrials 2]);
        obj.absorptions = [];
        absorptionsAllTrials = obj.compute(...
            oiSequence.frameAtIndex(idx), ...
            'emPath', emSubPath, ...
            'newNoise', newNoise, ...
            'currentFlag', false ...
            );
        
        % Next, compute full absorptions for the remaining OIs
        if (numel(indices)>1)
            for k = 2:numel(indices)
                % Update the @coneMosaic with the partial integration time
                if (k < numel(indices))
                    % full absorption during the duration of the OI
                    obj.integrationTime = oiRefreshInterval;
                else
                    % partial absorption during the last OI
                    idx = indices(end);
                    integrationTimeForLastPartialAbsorption = emEnd - oiTimeAxis(idx);
                    obj.integrationTime = integrationTimeForLastPartialAbsorption;
                end
                % Update the sum of partial integration times
                actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
                % Compute absorptions
                emSubPath = reshape(emPaths(1:nTrials, emIndex,:), [nTrials 2]);
                obj.absorptions = [];
                absorptionsAllTrials = absorptionsAllTrials + obj.compute(...
                    oiSequence.frameAtIndex(indices(k)), ...
                    'emPath', emSubPath, ...
                    'newNoise', newNoise, ...
                    'currentFlag', false ...
                    );
            end  % for k
        end
        
        % Check the sum of partial absorptions equal the default absorption time
        if (abs(actualIntegrationTime-defaultIntegrationTime) > defaultIntegrationTime*0.0001) && (emIndex < numel(eyeMovementTimeAxis))
            error('Actual integration time (%3.5f) not equal to desired value (%3.5f) [emIndex: %d / %d]\n', actualIntegrationTime, defaultIntegrationTime, emIndex, numel(eyeMovementTimeAxis));
        end
        
        % Reformat and insert to time series
        insertionIndices = round((eyeMovementTimeAxis(emIndex)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
        absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
    end % emIndex
end % oiRefreshInterval > defaultIntegrationTime

% Restore default integrationTime
obj.integrationTime = defaultIntegrationTime;

% Reload the full eye movement sequence for the first instance only
if (ndims(emPaths) == 3), obj.emPositions = squeeze(emPaths(1,:,:));
else                      obj.emPositions = emPaths;
end

% If we don't compute the current, we return only the absorptions from the
% last trial.  We never compute the current if there are no eye movements.
if (~currentFlag) || (numel(eyeMovementTimeAxis) == 1)
   
    if (isa(obj, 'coneMosaicHex'))
        % Return the absorptions from the last triale after reshaping to full 3D matrix [cone_rows, cone_cols, time]
        obj.absorptions = obj.reshapeHex2DmapToHex3Dmap(squeeze(absorptions(nTrials,:,:)));
    else
        % Reshape to full 4D matrix [instances, cone_rows, cone_cols, time]
        absorptions = reshape(absorptions, [nTrials size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
        
        % Return the absorptions from the last trial.
        obj.absorptions = reshape(squeeze(absorptions(nTrials,:,:)), [size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
    end
    return;
end

%% Photocurrent computation

% The currentFlag must be on, and there must be a few eye movements. So we
% compute. 
%

if (currentFlag)
    if (isa(obj, 'coneMosaicHex'))
        photocurrents = zeros(nTrials, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');
        for ii=1:nTrials
            % Reshape to full 3D matrix for obj.computeCurrent
            obj.absorptions = obj.reshapeHex2DmapToHex3Dmap(squeeze(absorptions(ii,:,:)));
            obj.computeCurrent;
            % Back to 2D matrix to save space
            photocurrents(ii,:,:) = single(obj.reshapeHex3DmapToHex2Dmap(obj.current));
        end
    else
        photocurrents = zeros(nTrials, obj.rows, obj.cols, numel(eyeMovementTimeAxis), 'single');
        for ii=1:nTrials
            obj.absorptions = reshape(squeeze(absorptions(ii,:,:,:)), [obj.rows obj.cols, numel(eyeMovementTimeAxis)]);
            obj.computeCurrent;
            photocurrents(ii,:,:,:) = single(reshape(obj.current, [1 obj.rows obj.cols numel(eyeMovementTimeAxis)]));
        end
    end
end


% Reload the absorptions from the last instance (again, since we destroyed
% obj.absorptions in the current computation)
if (isa(obj, 'coneMosaicHex'))
    % Return the absorptions from the last triale after reshaping to full 3D matrix [cone_rows, cone_cols, time]
    obj.absorptions = obj.reshapeHex2DmapToHex3Dmap(squeeze(absorptions(nTrials,:,:)));
else
    % Reshape to full 4D matrix [instances, cone_rows, cone_cols, time]
    absorptions = reshape(absorptions, [nTrials obj.rows obj.cols numel(eyeMovementTimeAxis)]);
    
    % Return the absorptions from the last trial.
    obj.absorptions = reshape(squeeze(absorptions(nTrials,:,:)), [obj.rows obj.cols numel(eyeMovementTimeAxis)]);
end

%% Nested function to reformat absorptions
    function reformatAbsorptionsAllTrialsMatrix(nTrials, timePointsNum, coneRows, coneCols)
        % Save all absorptions as singles to save space
        absorptionsAllTrials = single(absorptionsAllTrials);
        
        % Reshape to cones x instances x timePoints. Note the 3rd dimension of
        % absorptionsAllTrials is traditionally time, but when we are doing multiple instances, it is instances * time
        absorptionsAllTrials = reshape(absorptionsAllTrials, [coneRows*coneCols nTrials timePointsNum]);
        
        % Only get the absorptions for the non-null cones - this has an effect only for coneMosaicHex mosaics
        absorptionsAllTrials = absorptionsAllTrials(nonNullConesIndices,:,:);
        
        % Reshape to [instances x cones x timePoints]
        absorptionsAllTrials = permute(absorptionsAllTrials, [2 1 3]);
    end
end

%%
function displayProgress(workerID, workDescription, progress)

maxStarsNum = 60;
if (isnan(progress))
    fprintf('worker-%02d: %s |', workerID, workDescription);
    for k = 1:60
        fprintf('*');
    end
    fprintf('|');
else
    fprintf('worker-%02d: %s |', workerID, workDescription);
    if (progress>1)
        progress = 1;
    end
    starsNum = round(maxStarsNum*progress);
    for k = 1:starsNum
        fprintf('*');
    end
end
fprintf('\n');
end

