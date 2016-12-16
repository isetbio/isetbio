function [absorptions, photocurrents, LMSfilters] = computeForOISequence(obj, oiSequence, varargin)
% Compute cone absorptions and optionally photocurrents for a @oiSequence
%
%    [absorptions, photocurrents, LMSfilters] = cMosaic.compute(oiSequence, varargin);
%
% The oiSequence describes a temporal sequence of optical images.  This
% calculation converts the sequence into cone absorptions.  The algorithm
% is written to allow efficient calculation for multiple trials of the same
% sequence, but with different eye movement patterns.
%
% Inputs:
%   obj         - @coneMosaic object
%   oiSequence  - @oiSequence object
%
% Optional inputs:
%   emPaths      - [N x M x 2] matrix of N trials, each with a different
%                  Mx2 eye movement path
%   currentFlag  - logical, whether to compute photocurrent
%   seed         - noise seed
%
% Outputs:
%   absorptions          - cone photon absorptions (photon counts in integrationTime)
%   photocurrent         - cone photocurrent
%
% The simplest calculation is to send in a single oiSequence and a single
% eye movement sequence, stored in coneMosaic.emPositions
%
%   coneMosaic.compute(oiSequence)
%
% It is also possible to repeat the calculation for trials with different
% eye movement paths. First calculate the paths
%
%   tSamples = 100; nTrials = 10;
%   emPaths = coneMosaic.emGenSequence(tSamples,'nTrials',nTrials);
%
% Then send them in as an argument
%
%   absorptions, photocurrents] = coneMosaic.compute(oiSequence,'emPaths',emPaths)
%
% Upon return the coneMosaic object contains the eye movement positions
% from the last trial.  The absorption data for all the trials are
% contained in the returned outputs, absorptions and photocurrent, which
% are matrices that are [nTrials x row x col x time].
%
% You can control the photon noise by cmosaic.noiseFlag, and the
% photocurrent noise by cmosaic.os.noiseFlag.  These have the options
%
%    'random','frozen','none'
%
% We should allow that when 'frozen', you can send in a 'seed'.  May not be
% fully implemented yet.
%
% Examples:
%  This is an example of how to do this for 1,000 eye movement paths
%  (trials), each consisting of 100 eye movements.
%
%   nTrials = 50; nEyeMovements = 75; 
%   theEMPaths = cMosaic.emGenSequence(nEyeMovements,'nTrials',nTrials);
%   [absorptions, photocurrents] = cMosaic.computeForOISequence(...
%       theOIsequence, ...
%       'emPaths', theEMPaths, ...
%       'currentFlag', false);
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
% The same is true for the photocurrent return.
%
% See also: coneMosaic.compute, v_cmosaic, 
%
% NPC ISETBIO Team 2016

%% Parse inputs
p = inputParser;
p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
p.addParameter('seed',1, @isnumeric);                   % Seed for frozen noise
p.addParameter('emPaths', [], @isnumeric);
p.addParameter('currentFlag', false, @islogical);

p.addParameter('workerID', [], @isnumeric);
p.addParameter('workDescription', '', @ischar);
p.parse(oiSequence, varargin{:});

currentSeed     = p.Results.seed;
emPaths         = p.Results.emPaths;
currentFlag     = p.Results.currentFlag;
workerID        = p.Results.workerID;
workDescription = p.Results.workDescription;
oiTimeAxis      = oiSequence.timeAxis;

nTimes = numel(oiTimeAxis);
if (oiSequence.length ~= nTimes)
    error('oiTimeAxis and oiSequence must have equal length\n');
end

if (isempty(emPaths))
    % In this case we have an OI sequence, and we run with the one trial of
    % eye movement positions
    emPaths = reshape(obj.emPositions, [1 size(obj.emPositions)]);
end

if (isempty(emPaths))
    error('Either supply an ''emPaths'' key-value pair, or preload coneMosaic.emPositions');
end

%% Get ready for output variables

photocurrents = [];

%% Save default integration time
defaultIntegrationTime = obj.integrationTime;

%% Compute eye movement time axis
nTrials         = size(emPaths,1);
nEyeMovements = size(emPaths,2);
eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(nEyeMovements-1)) * obj.integrationTime;

%% Compute OIrefresh
if (numel(oiTimeAxis) == 1)
    % This isn't a sequence.  Why are we here at all?  Maybe we should
    % throw an error here.
    oiRefreshInterval = defaultIntegrationTime;
else
    % This is the time step of the oi sequence.
    oiRefreshInterval = oiSequence.timeStep; 
end


% Only allocate memory for the non-null cones in a 3D matrix 
% [instances x numel(nonNullConesIndices) x time].
% The LMS positions are indicated by pattern values of 2,3 or 4.  1 means
% blank (K).s
nonNullConesIndices = find(obj.pattern>1);
absorptions = zeros(nTrials, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');

% There are two main time sampling scenarios
if (oiRefreshInterval >= defaultIntegrationTime)
    % This one is when the oi update rate is SLOWER than the cone
    % integration time.  (Remember, that time is always equal to the eye
    % movement update rate.) 
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
    
    % We will make this optional in the future.  But BW wanted the waitbar
    % for help for a while.
    wb = waitbar(0,sprintf('Sequence %d',oiSequence.length));
    for oiIndex = 1:oiSequence.length
        % We pick out each image in the oiSequence.
        waitbar(oiIndex/oiSequence.length,wb);
        
        if (~isempty(workerID))
            % Update progress in command window if you have workerIDs. Not
            % sure what that is (BW).
            displayProgress(workerID, workDescription, 0.5*oiIndex/oiSequence.length);
        end
        
        % Current oi time limits
        tFrameStart = oiTimeAxis(oiIndex);
        tFrameEnd   = tFrameStart + oiRefreshInterval;
        
        % Find eye movement indices withing the oi limits
        indices = find( (eyeMovementTimeAxis >=  tFrameStart-defaultIntegrationTime) & ...
            (eyeMovementTimeAxis <= tFrameEnd-defaultIntegrationTime+eps) );
        
        if (isempty(indices))
            % Indices for the eye movement positions during this oi
            % presentation.
            error('Empty indices. This should never happen.');
        end
        
        % The first eye movement requires special treatment as it may have
        % started before the current frame, so we need to compute partial
        % absorptions over the previous frame and over the current frame
        idx = indices(1);
        integrationTimeForFirstPartialAbsorption  = tFrameStart - eyeMovementTimeAxis(idx);
        integrationTimeForSecondPartialAbsorption = eyeMovementTimeAxis(idx) + defaultIntegrationTime - tFrameStart;
        % fprintf('Integration times\n First %f\n Second %f\n',...
        %    integrationTimeForFirstPartialAbsorption,integrationTimeForSecondPartialAbsorption);
        % Partial absorptions (p1 in graph above) with previous oi
        % (across all instances)
        if (oiIndex > 1) % && integrationTimeForFirstPartialAbsorption > 1e-12
            % The usual case
            % Update the @coneMosaic with the partial integration time
            obj.integrationTime = integrationTimeForFirstPartialAbsorption;
            % Compute partial absorptions
            emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
            obj.absorptions = [];
            currentSeed = currentSeed  + 1;
            % Compute for all the eye movements, but just one frame
            absorptionsAllTrials = obj.compute(...
                oiSequence.frameAtIndex(oiIndex-1), ...
                'seed', currentSeed , ...
                'emPath', emSubPath, ...
                'currentFlag', false ...
                );
        else
            % We set to zero for the oiIndex == 1 condition
            absorptionsAllTrials = zeros(size(obj.pattern,1), size(obj.pattern,2), nTrials);
        end
        
        % Partial absorptions (p2 in graph above) with current oi across
        % all instances.
        % Update the @coneMosaic with the partial integration time
        if integrationTimeForSecondPartialAbsorption > 1e-12
            obj.integrationTime = integrationTimeForSecondPartialAbsorption;
            
            % Compute partial absorptions
            emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
            obj.absorptions = [];
            currentSeed = currentSeed  + 1;
            absorptionsDuringCurrentFrame =  obj.compute(...
                oiSequence.frameAtIndex(oiIndex), ...
                'seed', currentSeed, ...
                'emPath', emSubPath, ...
                'currentFlag', false ...
                );
            % vcNewGraphWin; imagesc(absorptionsDuringCurrentFrame);
            
            % summed absorptions
            absorptionsAllTrials = ...
                absorptionsAllTrials + absorptionsDuringCurrentFrame;
        end
        
        % Reformat and insert to time series
        insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
        absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
        
        % These are the indices of the next eye movements.
        % We compute the absorptions for this oi and default integration
        % time at the next eye movement positions.
        if (numel(indices)>1)
            % Update the @coneMosaic with the default integration time
            obj.integrationTime = defaultIntegrationTime;
            % Compute absorptions for all remaining OIs
            idx = indices(2:end);
            emSubPath = reshape(emPaths(1:nTrials, idx,:), [nTrials*numel(idx) 2]);
            obj.absorptions = [];
            currentSeed = currentSeed  + 1;
            absorptionsAllTrials = obj.compute(...
                oiSequence.frameAtIndex(oiIndex), ...
                'seed', currentSeed, ...
                'emPath', emSubPath, ...
                'currentFlag', false ...
                );
            
            % Reformat and insert to time series
            insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
            absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
        end
    end  % oiIndex
    
    delete(wb);

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
    for emIndex = 1:nEyeMovements
        
        if (~isempty(workerID))
            displayProgress(workerID, workDescription, 0.5*emIndex/nEyeMovements);
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
            disp('Fewer eye movement time samples than oi samples')
            break;
            %error('Empty indices. This should never happen.');
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
        currentSeed = currentSeed  + 1;
        absorptionsAllTrials = obj.compute(...
            oiSequence.frameAtIndex(idx), ...
            'seed', currentSeed, ...
            'emPath', emSubPath, ...
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
                currentSeed = currentSeed  + 1;
                absorptionsAllTrials = absorptionsAllTrials + ...
                    obj.compute(oiSequence.frameAtIndex(indices(k)), ...
                    'seed', currentSeed, ...
                    'emPath', emSubPath, ...
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
            currentSeed = currentSeed  + 1;
            LMSfilters = obj.computeCurrent('seed', currentSeed);
            % Back to 2D matrix to save space
            photocurrents(ii,:,:) = single(obj.reshapeHex3DmapToHex2Dmap(obj.current));
        end
    else
        photocurrents = zeros(nTrials, obj.rows, obj.cols, numel(eyeMovementTimeAxis), 'single');
        for ii=1:nTrials
            obj.absorptions = reshape(squeeze(absorptions(ii,:,:,:)), [obj.rows obj.cols, numel(eyeMovementTimeAxis)]);
            currentSeed = currentSeed  + 1;
            LMSfilters = obj.computeCurrent('seed', currentSeed);
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

